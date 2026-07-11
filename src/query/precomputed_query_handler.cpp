#include "query/precomputed_query_handler.h"
#include "common/utils.h"
#include <algorithm>

PrecomputedQueryHandler::PrecomputedQueryHandler(
    const NpyArray& edges_arr,
    const NpyArray& offsets_arr,
    const NpyArray& ts_arr,
    const NpyArray* normals_arr
){
    if (edges_arr.shape.size() != 2 || edges_arr.shape[1] != 2)
        throw std::runtime_error("PrecomputedQueryHandler: edges must be (E, 2)");
    if (offsets_arr.shape.size() != 1)
        throw std::runtime_error("PrecomputedQueryHandler: isect_offsets must be (E+1,)");

    size_t num_edges = edges_arr.shape[0];
    if (offsets_arr.shape[0] != num_edges + 1)
        throw std::runtime_error("PrecomputedQueryHandler: isect_offsets must have E+1 entries");

    const int32_t* edges = edges_arr.as_int32();
    const int32_t* offsets = offsets_arr.as_int32();
    const double* ts = ts_arr.as_float64();

    has_normals = (normals_arr != nullptr);
    const double* norms = nullptr;
    if (has_normals) {
        if (normals_arr->shape.size() != 2 || normals_arr->shape[1] != 3)
            throw std::runtime_error("PrecomputedQueryHandler: isect_normals must be (N, 3)");
        norms = normals_arr->as_float64();
    }

    edge_map.reserve(num_edges);
    for (size_t e = 0; e < num_edges; ++e) {
        size_t i = (size_t)edges[e * 2 + 0];
        size_t j = (size_t)edges[e * 2 + 1];
        if (i > j)
            throw std::runtime_error("PrecomputedQueryHandler: edges must be stored with smaller index first");

        int32_t start = offsets[e];
        int32_t end = offsets[e + 1];
        EdgeData data;
        data.ts.assign(ts + start, ts + end);
        if (has_normals) {
            data.normals.resize(end - start);
            for (int32_t k = start; k < end; ++k)
                data.normals[k - start] = Vector3{norms[k*3], norms[k*3+1], norms[k*3+2]};
        }
        edge_map[{i, j}] = std::move(data);
    }
}

void PrecomputedQueryHandler::query_intersections(
    const std::array<Vector3,4>& /*tet_positions*/,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool /*useRobust*/,
    bool /*recordNormals*/
){
    // Position-only overload cannot look up by index — return empty.
    // The pipeline should call the index-aware overload instead.
    for (int e = 0; e < 6; ++e) {
        edge_isect_ts[e].clear();
        edge_isect_normals[e].clear();
    }
}

void PrecomputedQueryHandler::query_intersections(
    const std::array<size_t,4>& tet_indices,
    const std::array<Vector3,4>& /*tet_positions*/,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool /*useRobust*/,
    bool /*recordNormals*/
){
    for (int e = 0; e < 6; ++e) {
        size_t local_i = ALL_TET_PAIRS[e].first;
        size_t local_j = ALL_TET_PAIRS[e].second;
        size_t global_i = tet_indices[local_i];
        size_t global_j = tet_indices[local_j];

        // Canonicalize: lookup key is (min, max)
        bool reversed = (global_i > global_j);
        size_t key_lo = reversed ? global_j : global_i;
        size_t key_hi = reversed ? global_i : global_j;

        auto it = edge_map.find({key_lo, key_hi});
        if (it == edge_map.end()) {
            edge_isect_ts[e].clear();
            edge_isect_normals[e].clear();
            continue;
        }

        const EdgeData& data = it->second;

        if (!reversed) {
            edge_isect_ts[e] = data.ts;
            edge_isect_normals[e] = data.normals;
        } else {
            // Reverse order and flip t-values: t -> 1-t
            size_t n = data.ts.size();
            edge_isect_ts[e].resize(n);
            for (size_t k = 0; k < n; ++k)
                edge_isect_ts[e][k] = 1.0 - data.ts[n - 1 - k];

            if (has_normals) {
                edge_isect_normals[e].resize(n);
                for (size_t k = 0; k < n; ++k)
                    edge_isect_normals[e][k] = data.normals[n - 1 - k];
            } else {
                edge_isect_normals[e].clear();
            }
        }
    }
}


void PrecomputedQueryHandler::query_normal(
    const Vector3& /*q*/, Vector3& /*normal*/, bool /*verbose*/
){
    throw std::runtime_error("PrecomputedQueryHandler: single-point normal query not supported");
}
