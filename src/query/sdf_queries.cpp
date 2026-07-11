#include "query/sdf_queries.h"


Vector3
grad_sdf(
    Vector3 q,
    const std::function<float(const Vector3&)>& sdf,
    float h
) {
    const Vector3 ex{h, 0, 0}, ey{0, h, 0}, ez{0, 0, h};
    const float fx1 = sdf(q + ex), fx0 = sdf(q - ex);
    const float fy1 = sdf(q + ey), fy0 = sdf(q - ey);
    const float fz1 = sdf(q + ez), fz0 = sdf(q - ez);
    return Vector3{(fx1 - fx0), (fy1 - fy0), (fz1 - fz0)} / (2.0f * h);
}


float linear_isect_t_interpolate(float t0, float t1, float sdf0, float sdf1){
    float abs_s0 = std::abs(sdf0);
    float abs_s1 = std::abs(sdf1);
    if (abs_s0 + abs_s1 < 1e-6f){
        return 0.5f * (t0 + t1);
    }
    return t0 + (t1 - t0) * (abs_s0 / (abs_s0 + abs_s1));
}

void edge_intersects_SDF(
    const Vector3& p0,
    const Vector3& p1,
    const std::function<float(const Vector3&)>& sdf,
    std::vector<float>& out_isect_ts,
    std::vector<Vector3>& out_isect_normals,
    size_t &query_count,
    float min_t_step
){
    out_isect_ts.clear();
    out_isect_normals.clear();
    query_count = 0;

    const size_t MAX_STEPS = 100000;

    float edge_length = (p1 - p0).norm();

    auto eval_point = [&](float t) -> Vector3 {
        return p0 + t * (p1 - p0);
    };

    auto sign_differs = [](float a, float b) -> bool {
        return (a > 0.0f && b < 0.0f) || (a < 0.0f && b > 0.0f);
    };

    float t = 0.0f;
    Vector3 p = eval_point(t);
    float sdf_val = sdf(p);

    // Track last non-zero sample so exact-zero landings don't trigger false sign changes
    float last_nonZ_sdf = sdf_val;
    float last_nonZ_t = (sdf_val != 0.0f) ? t : -1.0f;

    size_t step_count = 0;
    while (true){
        step_count++;
        if (step_count > MAX_STEPS){
            std::cerr << " SDF queries: Reached max steps along edge." << std::endl;
            return;
        }

        float next_t = t + std::max(std::abs(sdf_val) / edge_length, min_t_step);
        next_t = std::min(next_t, 1.0f);
        Vector3 next_p = eval_point(next_t);
        float next_sdf_val = sdf(next_p);
        query_count++;

        if (sign_differs(last_nonZ_sdf, next_sdf_val) && last_nonZ_t >= 0.0f){
            float isect_t = linear_isect_t_interpolate(last_nonZ_t, next_t, last_nonZ_sdf, next_sdf_val);
            Vector3 isect_p = eval_point(isect_t);
            Vector3 isect_normal = grad_sdf(isect_p, sdf).normalize();
            out_isect_ts.push_back(isect_t);
            out_isect_normals.push_back(isect_normal);
        }

        t = next_t;
        if (t >= 1.0f) break;
        sdf_val = next_sdf_val;
        if (sdf_val != 0.0f){
            last_nonZ_sdf = sdf_val;
            last_nonZ_t = t;
        }
    }
}


void tet_edge_intersections_SDF(
    const std::array<Vector3,4>& tet_positions,
    const std::function<float(const Vector3&)>& sdf,
    std::array<std::vector<double>,6>& out_edge_isect_ts,
    std::array<std::vector<Vector3>,6>& out_edge_isect_normals,
    std::array<size_t, 6> &query_count_per_edge,
    float tol
){
    for (size_t edge_idx = 0; edge_idx < 6; ++edge_idx){
        query_count_per_edge[edge_idx] = 0;
        int i = ALL_TET_PAIRS[edge_idx].first;
        int j = ALL_TET_PAIRS[edge_idx].second;
        std::vector<float> edge_isect_ts;
        size_t query_count_ij = 0;
        edge_intersects_SDF(
            tet_positions[i], tet_positions[j],
            sdf,
            edge_isect_ts, out_edge_isect_normals[edge_idx],
            query_count_ij,
            tol
        );
        query_count_per_edge[edge_idx] = query_count_ij;
        for (float t_val : edge_isect_ts) {
            out_edge_isect_ts[edge_idx].push_back(static_cast<double>(t_val));
        }
    }
}


// ============================================================================
// Wrapper functions for InputQueryHandler interface
// ============================================================================

Vector3 query_normal_SDF(
    const Vector3& q,
    const std::function<float(const Vector3&)>& sdf
) {
    return grad_sdf(q, sdf).normalize();
}
