
#include "query/mod2_reduction.h"

#include <numeric>


void convert_mod2_intersections(
    const std::array<std::vector<double>, 6>& in_edge_isects,
    std::array<std::vector<double>, 6>& out_mod2_isects,
    std::array<int, 6>& out_isect_counts
){
    for (size_t e = 0; e < 6; ++e){
        int count_isect = in_edge_isects[e].size();
        if (count_isect % 2 == 1){
            out_isect_counts[e] = 1;
            double avg_t = 0.0;
            for (double t : in_edge_isects[e])
                avg_t += t;
            avg_t /= static_cast<double>(count_isect);
            out_mod2_isects[e] = {avg_t};
        }
        else {
            out_isect_counts[e] = 0;
            out_mod2_isects[e] = {};
        }
    }
}


std::array<std::vector<Vector3>, 6>
convert_mod2_normals(
    const std::array<std::vector<double>, 6>& orig_ts,
    const std::array<std::vector<double>, 6>& mod2_ts,
    const std::array<std::vector<Vector3>, 6>& orig_normals
){
    std::array<std::vector<Vector3>, 6> result;
    for (size_t e = 0; e < 6; ++e) {
        if (mod2_ts[e].empty()) continue;
        double target = mod2_ts[e][0];
        size_t best = 0;
        double best_dist = std::abs(orig_ts[e][0] - target);
        for (size_t k = 1; k < orig_ts[e].size(); ++k) {
            double d = std::abs(orig_ts[e][k] - target);
            if (d < best_dist) { best_dist = d; best = k; }
        }
        result[e] = {orig_normals[e][best]};
    }
    return result;
}


bool apply_mod2_reduction(
    std::array<std::vector<double>,6>& isect_ts,
    std::array<int,6>& isect_counts
){
    std::array<std::vector<double>,6> mod2_ts;
    convert_mod2_intersections(isect_ts, mod2_ts, isect_counts);
    if (std::accumulate(isect_counts.begin(), isect_counts.end(), 0) == 0)
        return true;
    isect_ts = std::move(mod2_ts);
    return false;
}


bool apply_mod2_reduction(
    std::array<std::vector<double>,6>& isect_ts,
    std::array<std::vector<Vector3>,6>& isect_normals,
    std::array<int,6>& isect_counts
){
    std::array<std::vector<double>,6> mod2_ts;
    convert_mod2_intersections(isect_ts, mod2_ts, isect_counts);
    if (std::accumulate(isect_counts.begin(), isect_counts.end(), 0) == 0)
        return true;
    isect_normals = convert_mod2_normals(isect_ts, mod2_ts, isect_normals);
    isect_ts = std::move(mod2_ts);
    return false;
}
