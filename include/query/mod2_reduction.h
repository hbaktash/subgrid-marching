#pragma once

#include "common/utils.h"

// Reduce intersections to mod-2: odd-count edges get a single averaged t,
// even-count edges get cleared. Updates out_isect_counts in place.
void convert_mod2_intersections(
    const std::array<std::vector<double>, 6>& edge_intersections,
    std::array<std::vector<double>, 6>& out_mod2_intersections,
    std::array<int, 6>& out_isect_counts
);

// For each mod2-reduced edge, pick the original normal whose t is closest
// to the mod2 average. Even-count edges (empty in mod2_ts) get empty output.
std::array<std::vector<Vector3>, 6>
convert_mod2_normals(
    const std::array<std::vector<double>, 6>& orig_ts,
    const std::array<std::vector<double>, 6>& mod2_ts,
    const std::array<std::vector<Vector3>, 6>& orig_normals
);

// Apply mod2 reduction to t-values (primal pipeline, no normals).
// Returns true if the tet should be skipped (all edges reduced to zero).
bool apply_mod2_reduction(
    std::array<std::vector<double>,6>& isect_ts,
    std::array<int,6>& isect_counts
);

// Apply mod2 reduction to both t-values and normals (dual pipeline).
// Returns true if the tet should be skipped (all edges reduced to zero).
bool apply_mod2_reduction(
    std::array<std::vector<double>,6>& isect_ts,
    std::array<std::vector<Vector3>,6>& isect_normals,
    std::array<int,6>& isect_counts
);
