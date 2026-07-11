#pragma once

#include "nc/NC_solver.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/cv_interpolation.h"


using std::vector;
using std::array;
using std::tuple;
using std::pair;


// helpers (get_interior_vertex and extract_curve_edge_isect_ts are static in the .cpp)

EdgeInts
combface_to_edgeints(
    const CombFace &comb_face
);

vector<int>
separated_vertices(
    const EdgeInts &cruve_edge_ints
);

// construction
TriangleSoup
evensum_surface_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &global_edge_isect_ts,
    const vector<CombFace> &non_normal_comb_curves,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
);


TriangleSoup
diagonal_surface_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &global_edge_isect_ts,
    const CombFace &diagonal_comb_curve,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
);


TriangleSoup
nondiagonal_surface_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &global_edge_isect_ts,
    const CombFace &nondiagonal_comb_curve,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
);


// construction stages, called in order by nondiagonal_surface_construction.
// these are combinatorial only: they build CombFaces; the caller adds scoop
// mid vertices and emits them into the soup.
vector<CombFace>
smeared_quads_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex
);

vector<CombFace>
interior_region_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex
);

vector<CombFace>
corner_triangle_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex,
    bool scoop_mid_vertices
);

vector<CombFace>
tunnel_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex
);
