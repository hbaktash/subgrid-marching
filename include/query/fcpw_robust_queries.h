#pragma once

#include "fcpw/fcpw.h"

#include <vector>

#ifdef HAVE_FCPW

// Repeatedly applies FCPW's robust single-hit ray query to collect all hits along a ray.
// The ray is advanced by a small offset after each hit to avoid returning the same
// intersection repeatedly at shared boundaries or numerically unstable contacts.
std::vector<fcpw::Interaction<3>> intersect_robust_all_hits(
    fcpw::Scene<3>& accel,
    fcpw::Ray<3> ray,
    bool recordNormals = true,
    float advanceEpsilon = 1e-5f,
    int maxHits = -1
);

#endif