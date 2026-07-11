#pragma once

#include "nc/NC_solver.h"

using std::vector;
using std::array;
using std::tuple;
using std::pair;


// returns lists of cobinatorial boundary polygons; no positions
// - 3 lists of CombinatorialFaces; each CombFace is a list of CombinatorialVertex a type
//    - a: open curves
//    - b: non-normal closed curves
//    - c: normal closed curves
// optionally returns the tet face that each segment sits on for Dual use
tuple<vector<CombFace>, vector<CombFace>, vector<CombFace>>
boundary_comb_curves(
    const array<size_t,4> tet_global_indices,
    const array<int,6> edge_isect_counts,
    EdgeOccupations &out_edge_occupations,
    bool populate_normal_curves
);


bool pick_residual_cv(
    int i, int j, int other,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds,
    CombinatorialVertex &out_cv
);

void pick_residuals_and_scoop_cvs(
    const NormalCoordinates &NC,
    const array<size_t,4> &tet_global_indices,
    vector<CombinatorialVertex> &residual_cvs,
    vector<CombinatorialVertex> &scoop_cvs
);

// pickout the odd vertex and scoop vertices on edge {i,j} 
std::vector<CombinatorialVertex>
pick_scoop_cvs(
    int i, int j, int other,
    const NormalCoordinates &nc,
    const CombinatorialVertex &residual_cv
);


void 
trace_generic_curve(
    const CombinatorialVertex &start_cv,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds,
    EdgeOccupations &occupations,
    std::vector<CombinatorialVertex> &curve_vertices,
    bool keep_occupation_track = true
);



CombinatorialVertex
find_next_generic_vertex(
    const CombinatorialVertex &cv0,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds
);


// Returns the local tet face {edge_i, edge_j, third_vertex} for segment v_a → v_b.
// Matches the {cv_next.i, cv_next.j, cv_next.k} convention from the dual pipeline.
std::array<size_t,3>
segment_tet_face(
    const CombinatorialVertex &v_a,
    const CombinatorialVertex &v_b
);

// For each segment in comb_face, returns the local tet face {i,j,k}.
// Closed curves: n vertices → n segments (last wraps v_{n-1} → v_0).
// Open curves:   n vertices → n-1 segments.
// Output index i corresponds to segment (vertices[i] → vertices[(i+1) % n]).
std::vector<std::array<size_t,3>>
curve_segment_faces(
    const CombFace &comb_face
);
