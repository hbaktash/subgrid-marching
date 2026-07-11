#pragma once

#include "nc/NC_solver.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/cv_interpolation.h"
#include "subgrid_MT/normal_reconstruction.h"
#include "subgrid_MT/evensum_construction.h"

using std::vector;
using std::array;
using std::tuple;
using std::pair;


std::array<std::vector<double>,6> 
unoccupied_edge_isect_ts(
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const EdgeOccupations &edge_occupations
);


TriangleSoup
subgrid_surface(
    const std::array<Vector3,4> &tet_positions,
    const std::array<size_t,4> &tet_global_indices,
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const vector<CombFace> &open_boundary_curves,
    const vector<CombFace> &scoop_boundary_curves,
    const EdgeOccupations &non_normal_edge_occupations,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge,
    bool verbose = false
);

TriangleSoup
subgrid_greedy(
    const std::array<Vector3,4> &tet_positions,
    const std::array<size_t,4> &tet_global_indices,
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const vector<CombFace> &scoop_boundary_curves,
    const vector<CombFace> &normal_boundary_curves,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
);

