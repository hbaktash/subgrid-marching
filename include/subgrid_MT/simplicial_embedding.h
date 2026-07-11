#pragma once

#include "common/utils.h"
#include "nc/NC_solver.h"


using std::vector;
using std::array;


// for two cv's that are on the same edge, and are neighbors;
// make sure i and j are the same, and edge_int is the same; order can be different by 1
void
align_cvs(
    const CombinatorialVertex &cv1,
    const CombinatorialVertex &cv2,
    CombinatorialVertex &aligned_cv1,
    CombinatorialVertex &aligned_cv2
);

// check if the inner segment between two combinatorial vertices
// is capped by a scoop on either side, and return capped sides
vector<int>
get_close_scoop_sides(
    const NormalCoordinates &smeared_curve_nc,
    const CombinatorialVertex &inner_segment_cv1,
    const CombinatorialVertex &inner_segment_cv2
);


// SPIRAL faces only: walks the curve's own vertices, no NormalCoordinates needed.
// Call this directly when the face is known to be a spiral (e.g. diagonal curves).
CombFace
add_mid_scoop_vertices_spiral(
    const CombFace &spiral_comb_face
);

// adds mid scoop steiner vertices to an even-sum comb face, dispatching on its type
CombFace
add_mid_scoop_vertices(
    const CombFace &even_sum_comb_face,
    const NormalCoordinates &smeared_curve_nc
);
