#pragma once

#include "common/triangle_soup.h"
#include "subgrid_MT/boundary_curve.h"

// Per-tet dual subgrid construction: from boundary curves + intersection data,
// builds a TriangleSoup with positions, faces, faces_per_edge, signatures, and
// dual_positions (QEF-solved), remapped to global tet indices.
TriangleSoup
dual_subgrid_surface(
    const std::array<Vector3,4>& tet_positions,
    const std::array<size_t,4>& tet_indices,
    const std::array<std::vector<double>,6>& edge_isect_ts,
    const std::array<std::vector<Vector3>,6>& edge_isect_normals,
    const std::vector<CombFace>& open_curves,
    const std::vector<CombFace>& scoop_curves,
    const std::vector<CombFace>& normal_curves,
    double reg_alpha,
    bool project_duals
);


// Extract normals (parallel to CVs across all curves) from edge_isect_normals.
std::vector<Vector3>
normals_from_cv_polygons(
    const std::vector<CombFace>& curves,
    const std::array<std::vector<double>,6>& edge_isect_ts,
    const std::array<std::vector<Vector3>,6>& edge_isect_normals
);

// Build a local TriangleSoup from boundary curves: interpolates positions,
// populates NORMAL-type signatures and faces_per_edge.
TriangleSoup
build_dual_local_soup(
    const std::vector<CombFace>& curves,
    const std::array<Vector3,4>& tet_positions,
    const std::array<std::vector<double>,6>& edge_isect_ts
);
