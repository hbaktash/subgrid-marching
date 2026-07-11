#pragma once

#include "nc/NC_solver.h"
#include "common/triangle_soup.h"

using std::vector;
using std::array;
using std::pair;


// ---- Curve stuff; mainly for visuals ----
vector<Vector3>
comb_curve_to_vec3_curve(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    double bulge
);


// ---- Geometry embedding ----
Vector3
interpolate_comb_vertex(
    const CombinatorialVertex &cv,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts,
    double bulge
);

Vector3
interpolate_derived_vertex(
    const DerivedVertex &dv,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts
);

TriangleSoup
embed_comb_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    bool mid_scoop_vertices,
    double bulge
);
