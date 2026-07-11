#pragma once

#include "nc/NC_solver.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/cv_interpolation.h"


using std::vector;
using std::array;
using std::tuple;
using std::pair;


TriangleSoup
normal_surface_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    bool verbose
);


// Normal stuff
TriangleSoup
corner_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
);


TriangleSoup
diagonal_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
);


// Almost Normal stuff
TriangleSoup
octagon_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
);


TriangleSoup spiral_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
);


TriangleSoup recursive_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC,
    bool verbose
);


std::pair<Vector3, std::array<double, 4>>
compute_recursion_center_and_isect_ts(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
);