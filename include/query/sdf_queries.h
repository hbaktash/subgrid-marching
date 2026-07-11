#pragma once

#include "common/utils.h"

Vector3
grad_sdf(
    Vector3 q,
    const std::function<float(const Vector3&)>& sdf,
    float h = 1e-5
);

void edge_intersects_SDF(
    const Vector3& p0,
    const Vector3& p1,
    const std::function<float(const Vector3&)>& sdf,
    std::vector<float>& out_isect_ts,
    std::vector<Vector3>& out_isect_normals,
    size_t &query_count,
    float min_step
);

void tet_edge_intersections_SDF(
    const std::array<Vector3,4>& tet_positions,
    const std::function<float(const Vector3&)>& sdf,
    std::array<std::vector<double>,6>& out_edge_isect_ts,
    std::array<std::vector<Vector3>,6>& out_edge_isect_normals,
    std::array<size_t, 6> &query_count_per_edge,
    float min_step = 1e-3
);

Vector3 query_normal_SDF(
    const Vector3& q,
    const std::function<float(const Vector3&)>& sdf
);
