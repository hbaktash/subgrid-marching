#pragma once

#include "common/utils.h"

#include "geometrycentral/surface/meshio.h"


void get_stem_parent_extension(
    const std::string &full_filepath,
    std::string &stem,
    std::string &parent,
    std::string &extension
);

void 
ensure_path_exists(
    const std::string &full_filepath
);

void save_polygon_soup_as_obj(
    const std::string& filename,
    const std::vector<Vector3>& vertices,
    const std::vector<std::vector<size_t>>& polygons
);

struct TetConfigFileData {
    std::array<int,6> edge_ints;
    std::array<std::vector<double>,6> edge_isect_ts;
    std::vector<std::pair<int,int>> intersecting_faces; // optional, empty if not present
    double bulge = -1; // -1 means not specified in file
};

bool load_tet_config_file(const std::string& path, TetConfigFileData& out);