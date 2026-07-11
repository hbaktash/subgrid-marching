#include "query/fcpw_robust_queries.h"

#ifdef HAVE_FCPW

#include <algorithm>
#include <cmath>
#include <iostream>

std::vector<fcpw::Interaction<3>> intersect_robust_all_hits(
    fcpw::Scene<3>& accel,
    fcpw::Ray<3> ray,
    bool recordNormals,
    float advanceEpsilon,
    int maxHits
){
    std::vector<fcpw::Interaction<3>> intersections;

    if (advanceEpsilon <= 0.0f) {
        advanceEpsilon = 1e-5f;
    }

    const float initialTMax = ray.tMax;
    float traveled = 0.0f;
    int iteration = 0;

    while (ray.tMax > 0.0f) {
        float remainingDistance = ray.tMax;

        fcpw::Interaction<3> intersection;
        if (!accel.intersectRobust(ray, intersection)) {
            break;
        }
        if (!std::isfinite(intersection.d) || intersection.d < 0.0f) {
            break;
        }

        // Convert local ray distance into distance along the original segment.
        float globalDistance = traveled + intersection.d;
        if (!std::isfinite(globalDistance) || globalDistance < 0.0f || globalDistance > initialTMax + advanceEpsilon) {
            break;
        }

        intersection.d = globalDistance;

        if (!recordNormals) {
            intersection.n = fcpw::Vector3::Zero();
        }

        intersections.push_back(intersection);

        if (maxHits > 0 && static_cast<int>(intersections.size()) >= maxHits) {
            break;
        }

        float step = std::max(advanceEpsilon, globalDistance - traveled + advanceEpsilon);
        if (!std::isfinite(step) || step <= 0.0f || step >= remainingDistance) {
            break;
        }

        fcpw::Vector3 nextOrigin = ray.o + ray.d * step;
        traveled += step;
        ray = fcpw::Ray<3>(nextOrigin, ray.d, remainingDistance - step);
    }

    return intersections;
}

#endif
