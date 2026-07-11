#include "common/npz_reader.h"

#include <fstream>
#include <cstring>
#include <algorithm>

// NPZ is a ZIP archive of .npy files. We parse the ZIP local file headers
// directly (no deflation — NumPy always stores arrays uncompressed).
// NPY format: magic "\x93NUMPY" + version(2B) + header_len + python dict header + raw data.

size_t NpyArray::num_elements() const {
    size_t n = 1;
    for (auto s : shape) n *= s;
    return n;
}

const double* NpyArray::as_float64() const {
    if (!is_float || word_size != 8)
        throw std::runtime_error("NpyArray: requested float64 but array is not float64");
    return reinterpret_cast<const double*>(data.data());
}

const int32_t* NpyArray::as_int32() const {
    if (is_float || word_size != 4)
        throw std::runtime_error("NpyArray: requested int32 but array is not int32");
    return reinterpret_cast<const int32_t*>(data.data());
}

bool NpzFile::has(const std::string& name) const {
    return arrays.count(name) > 0;
}

const NpyArray& NpzFile::operator[](const std::string& name) const {
    auto it = arrays.find(name);
    if (it == arrays.end())
        throw std::runtime_error("NpzFile: array '" + name + "' not found");
    return it->second;
}

// Parse the dtype descriptor from the NPY header dict
static void parse_dtype(const std::string& descr, size_t& word_size, bool& is_float) {
    // descr looks like "'<f8'", "'<i4'", "'>f8'", etc.
    // We strip quotes and look at the type char and size
    std::string d = descr;
    d.erase(std::remove(d.begin(), d.end(), '\''), d.end());
    d.erase(std::remove(d.begin(), d.end(), ' '), d.end());

    // Skip endian character (< or > or |)
    size_t pos = 0;
    if (d[0] == '<' || d[0] == '>' || d[0] == '|') pos = 1;

    char type_char = d[pos];
    word_size = std::stoul(d.substr(pos + 1));

    is_float = (type_char == 'f');
    if (type_char != 'f' && type_char != 'i' && type_char != 'u')
        throw std::runtime_error("npz_reader: unsupported dtype '" + descr + "'");
}

// Parse shape tuple from the NPY header dict
static std::vector<size_t> parse_shape(const std::string& shape_str) {
    // shape_str looks like "(128, 3)" or "(5,)" or "(100,3,)"
    std::vector<size_t> shape;
    std::string s = shape_str;
    s.erase(std::remove(s.begin(), s.end(), '('), s.end());
    s.erase(std::remove(s.begin(), s.end(), ')'), s.end());
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());

    size_t start = 0;
    while (start < s.size()) {
        size_t comma = s.find(',', start);
        if (comma == std::string::npos) comma = s.size();
        std::string token = s.substr(start, comma - start);
        if (!token.empty())
            shape.push_back(std::stoul(token));
        start = comma + 1;
    }
    return shape;
}

// Extract a value for a given key from the numpy header dict string
static std::string extract_dict_value(const std::string& header, const std::string& key) {
    std::string search = "'" + key + "'";
    size_t pos = header.find(search);
    if (pos == std::string::npos)
        throw std::runtime_error("npz_reader: key '" + key + "' not found in header");
    pos = header.find(':', pos);
    if (pos == std::string::npos)
        throw std::runtime_error("npz_reader: malformed header");
    pos++;
    while (pos < header.size() && header[pos] == ' ') pos++;

    // Value can be a string (quoted), tuple (parenthesized), or identifier
    if (header[pos] == '\'') {
        size_t end = header.find('\'', pos + 1);
        return header.substr(pos, end - pos + 1);
    } else if (header[pos] == '(') {
        size_t end = header.find(')', pos);
        return header.substr(pos, end - pos + 1);
    } else {
        size_t end = header.find_first_of(",}", pos);
        return header.substr(pos, end - pos);
    }
}

static NpyArray parse_npy_from_buffer(const char* buf, size_t buf_size) {
    // Verify magic
    if (buf_size < 10 || std::memcmp(buf, "\x93NUMPY", 6) != 0)
        throw std::runtime_error("npz_reader: invalid NPY magic");

    uint8_t major = buf[6];
    // uint8_t minor = buf[7];

    uint32_t header_len;
    size_t header_offset;
    if (major == 1) {
        header_len = *reinterpret_cast<const uint16_t*>(buf + 8);
        header_offset = 10;
    } else {
        header_len = *reinterpret_cast<const uint32_t*>(buf + 8);
        header_offset = 12;
    }

    std::string header(buf + header_offset, header_len);
    size_t data_offset = header_offset + header_len;

    NpyArray arr;
    std::string descr = extract_dict_value(header, "descr");
    parse_dtype(descr, arr.word_size, arr.is_float);
    arr.shape = parse_shape(extract_dict_value(header, "shape"));

    size_t data_size = arr.num_elements() * arr.word_size;
    if (data_offset + data_size > buf_size)
        throw std::runtime_error("npz_reader: NPY data exceeds buffer");

    arr.data.resize(data_size);
    std::memcpy(arr.data.data(), buf + data_offset, data_size);
    return arr;
}

NpzFile load_npz(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file)
        throw std::runtime_error("npz_reader: cannot open '" + path + "'");

    NpzFile result;

    // Read ZIP local file headers
    while (true) {
        // ZIP local file header signature: PK\x03\x04
        char sig[4];
        if (!file.read(sig, 4)) break;
        if (std::memcmp(sig, "PK\x03\x04", 4) != 0) break;

        // Skip to filename length and extra length fields
        file.seekg(22, std::ios::cur); // skip rest of fixed header (26 bytes total - 4 sig = 22)
        uint16_t fname_len, extra_len;
        file.read(reinterpret_cast<char*>(&fname_len), 2);
        file.read(reinterpret_cast<char*>(&extra_len), 2);

        // Read filename
        std::string fname(fname_len, '\0');
        file.read(&fname[0], fname_len);

        // Skip extra field
        file.seekg(extra_len, std::ios::cur);

        // We need the compressed size to know how much data to read.
        // Go back and read it from the fixed header.
        auto current_pos = file.tellg();
        file.seekg(-(std::streamoff)(fname_len + extra_len + 4 + 22), std::ios::cur);
        // Now at offset 4 from start of this header. Skip to compressed_size at offset 18.
        file.seekg(14, std::ios::cur);
        uint32_t compressed_size;
        file.read(reinterpret_cast<char*>(&compressed_size), 4);
        // Return to data start
        file.seekg(current_pos);

        // Read the data (stored uncompressed for .npy in npz)
        std::vector<char> buf(compressed_size);
        file.read(buf.data(), compressed_size);

        // Strip .npy extension from filename
        std::string array_name = fname;
        if (array_name.size() > 4 && array_name.substr(array_name.size() - 4) == ".npy")
            array_name = array_name.substr(0, array_name.size() - 4);

        result.arrays[array_name] = parse_npy_from_buffer(buf.data(), buf.size());
    }

    if (result.arrays.empty())
        throw std::runtime_error("npz_reader: no arrays found in '" + path + "'");

    return result;
}
