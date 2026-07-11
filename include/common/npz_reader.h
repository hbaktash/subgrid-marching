#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <stdexcept>

// Minimal NPZ reader: loads named arrays from a .npz file.
// Supports float64 and int32 arrays only (sufficient for our explicit input format).

struct NpyArray {
    std::vector<char> data;
    std::vector<size_t> shape;
    size_t word_size = 0;   // bytes per element (4 for int32, 8 for float64)
    bool is_float = false;  // true = float64, false = int32

    size_t num_elements() const;

    // Typed accessors (throw if type mismatch)
    const double* as_float64() const;
    const int32_t* as_int32() const;
};

struct NpzFile {
    std::unordered_map<std::string, NpyArray> arrays;

    bool has(const std::string& name) const;
    const NpyArray& operator[](const std::string& name) const;
};

// Load all arrays from an .npz file
NpzFile load_npz(const std::string& path);
