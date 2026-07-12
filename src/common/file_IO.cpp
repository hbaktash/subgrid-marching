#include "common/file_IO.h"
#include <filesystem>
#include <fstream>
#include <sstream>


void get_stem_parent_extension(
    const std::string &full_filepath,
    std::string &stem,
    std::string &parent,
    std::string &extension
){
    std::filesystem::path filepath(full_filepath);
    stem = filepath.stem().string();
    parent = filepath.parent_path().string();
    extension = filepath.extension().string();
}


// helper for paths; if path does not exist, create it
void ensure_path_exists(const std::string &full_filepath) {
    std::filesystem::path filepath(full_filepath);
    std::filesystem::path dir = filepath.parent_path();
    if (!dir.empty() && !std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
}

bool load_tet_config_file(const std::string& path, TetConfigFileData& out) {
    std::ifstream f(path);
    if (!f.is_open()) return false;

    out.intersecting_faces.clear();
    std::string line;

    while (std::getline(f, line)) {
        if (line.find("edge_ints:") != std::string::npos) {
            std::istringstream ss(line.substr(line.find(':') + 1));
            for (int i = 0; i < 6; ++i) ss >> out.edge_ints[i];
        }
        else if (line.find("bulge:") != std::string::npos && line.find("bulge") < 6) {
            std::istringstream ss(line.substr(line.find(':') + 1));
            ss >> out.bulge;
        }
        else if (line.find("intersecting face pairs") != std::string::npos) {
            while (std::getline(f, line) && !line.empty() && line[0] == ' ') {
                std::istringstream ss(line);
                int a, b;
                if (ss >> a >> b) out.intersecting_faces.push_back({a, b});
            }
            // The line that broke the loop might be the t-values header
            if (line.find("t-values per edge") != std::string::npos) goto parse_tvalues;
        }
        else if (line.find("t-values per edge") != std::string::npos) {
            parse_tvalues:
            for (int ei = 0; ei < 6; ++ei) {
                if (!std::getline(f, line)) return false;
                size_t colon = line.find(':');
                if (colon == std::string::npos) return false;
                std::istringstream ss(line.substr(colon + 1));
                out.edge_isect_ts[ei].clear();
                double t;
                while (ss >> t) out.edge_isect_ts[ei].push_back(t);
            }
        }
    }
    return true;
}

void save_polygon_soup_as_obj(
    const std::string& filename,
    const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& polygons
){
    ensure_path_exists(filename);
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write vertices
    for (const auto& v : vertices) {
        out_file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    // Write faces
    for (const auto& poly : polygons) {
        out_file << "f";
        for (const auto& idx : poly) {
            out_file << " " << (idx + 1); // OBJ format is 1-indexed
        }
        out_file << "\n";
    }
    out_file.close();
}