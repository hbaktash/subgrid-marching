#pragma once 

#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

// polyscope for Debugging (only when built with the viewer)
#ifdef HAVE_POLYSCOPE
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#endif

#include <vector>
#include <set>
#include <initializer_list>
#include <unordered_set>
#include <iostream>
#include <string>


using namespace geometrycentral;
using namespace geometrycentral::surface;

using std::vector;

#define ANSI_FG_MAGENTA "\x1b[35m"
#define ANSI_FG_YELLOW "\x1b[33m"
#define ANSI_FG_GREEN "\x1b[32m"
#define ANSI_FG_WHITE "\x1b[37m"
#define ANSI_FG_RED "\x1b[31m"
#define ANSI_RESET "\x1b[0m"


// Lightweight tagged diagnostics with a consistent prefix + color.
// Only the tag is colored so the message text stays copy-pasteable.
// Info goes to stdout; warnings and errors go to stderr.
inline void log_info(const std::string& msg)  { std::cout << ANSI_FG_GREEN  << "[info] "  << ANSI_RESET << msg << "\n"; }
inline void log_warn(const std::string& msg)  { std::cerr << ANSI_FG_YELLOW << "[warn] "  << ANSI_RESET << msg << "\n"; }
inline void log_error(const std::string& msg) { std::cerr << ANSI_FG_RED    << "[error] " << ANSI_RESET << msg << "\n"; }



void print_progress(double progress);


void str_to_pair(const std::string& str, size_t& a, size_t& b);


bool triangle_inequality(int a, int b, int c, 
                         int &res_a, int &res_b, int &res_c);


//  set operations on {0,1,2,3}
std::vector<int> complement(std::initializer_list<int> inputs);
// operate on vector input too
std::vector<int> complement(const std::vector<int> &inputs);

std::vector<std::pair<int, int>> all_pairs(std::vector<int> universe);

std::vector<std::pair<int, int>> all_pairs(std::array<int, 3> universe);


inline std::array<std::pair<int, int>,6>
ALL_TET_PAIRS = {
    std::pair<int,int>{0,1}, std::pair<int,int>{0,2}, std::pair<int,int>{0,3},
    std::pair<int,int>{1,2}, std::pair<int,int>{1,3}, std::pair<int,int>{2,3}
};

// quick function to map edge pair to index
inline size_t edge_pair_to_index(int i, int j) {
    if (i > j) std::swap(i, j);
    if (i == 0 && j == 1) return 0;
    if (i == 0 && j == 2) return 1;
    if (i == 0 && j == 3) return 2;
    if (i == 1 && j == 2) return 3;
    if (i == 1 && j == 3) return 4;
    if (i == 2 && j == 3) return 5;
    std::cout << "Invalid edge pair: (" << i << ", " << j << ")" << std::endl;
    throw std::invalid_argument("edge_pair_to_index: Invalid edge pair");
};

inline std::array<std::array<int,3>,4>
ALL_TET_TRIPLETS = {
    std::array<int,3>{1,2,3},
    std::array<int,3>{0,2,3},
    std::array<int,3>{0,1,3},
    std::array<int,3>{0,1,2}
};

struct Tet{
    std::array<size_t,4> global_ind;
    // std::vector<std::vector<size_t>> faces; // trivial set of faces; orientation does not matter if global to local is assigned properly?
    std::array<Vector3,4> pos; // positions
};

void center_and_normalize(SurfaceMesh &mesh, VertexPositionGeometry &geometry, double r = 0.99);

void center_and_normalize(std::vector<Vector3> &points, double r = 0.99);

void triangulate_polygons(
    const std::vector<std::vector<size_t>> &input_polygons,
    std::vector<std::vector<size_t>> &output_triangles,
    bool remove_segments = true
);

std::vector<std::vector<size_t>> 
segments_to_polygons(
    const std::vector<std::pair<size_t, size_t>> &segments,
    const size_t pool_size
);

std::vector<size_t> 
pairs_to_single_polygon(
    const std::vector<std::pair<size_t, size_t>> &segments,
    bool warn = false
);

bool is_inside_tet(
    const Vector3 &p,
    const std::array<Vector3,4> &tet_positions,
    double eps = 1e-12
);


void
pairss_to_tupless(
    const std::vector<std::vector<std::pair<size_t, size_t>>> &pairss,
    std::vector<std::vector<std::tuple<size_t, size_t>>> &tupless
);


struct VoxelActivityStats {
    size_t active_voxels = 0;
    size_t active_or_neighbor_voxels = 0;
    size_t active_vertices = 0;
    size_t active_or_neighbor_vertices = 0;
    size_t active_voxel_equivalent_tets = 0;
};


class VoxelActivityTracker {
public:
    explicit VoxelActivityTracker(size_t grid_resolution);

    // Report tet activity using explicit tet vertex indices.
    void report_tet(const std::array<size_t, 4>& tet_indices, bool has_polygon);

    // Return current stats.
    VoxelActivityStats stats() const;

private:
    size_t grid_resolution_;

    std::unordered_set<size_t> active_voxels_;
    std::unordered_set<size_t> active_or_neighbor_voxels_;
    std::unordered_set<size_t> active_vertices_;
    std::unordered_set<size_t> active_or_neighbor_vertices_;

    size_t tet_indices_to_voxel_index(const std::array<size_t, 4>& tet_indices) const;
    std::array<size_t, 3> voxel_index_to_coords(size_t voxel_index) const;
    size_t voxel_coords_to_index(size_t x, size_t y, size_t z) const;
    size_t vertex_coords_to_index(size_t x, size_t y, size_t z) const;
    void mark_voxel_vertices(size_t voxel_index, std::unordered_set<size_t>& vertex_set);
    void mark_active_voxel_neighbors(size_t voxel_index);
};


bool segment_triangle_intersection(
    const geometrycentral::Vector3& p0,
    const geometrycentral::Vector3& p1,
    const geometrycentral::Vector3& a,
    const geometrycentral::Vector3& b,
    const geometrycentral::Vector3& c,
    double& out_t
);

// edge convention: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
bool check_even_sum(
    const std::array<int, 6> e_isects
);

bool check_triangle_ineq(
    const std::array<int, 6> e_isects
);

vector<double>
generate_uniform_isect_ts_single_edge(
    int edge_int, 
    double t_max = 1.0
);