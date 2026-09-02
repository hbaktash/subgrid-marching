#include "query/edge_isect_cache.h"
#include "common/parallel_for.h"
#include "common/utils.h"

#include <stdexcept>


void emit_edge_isect(
    const EdgeIsect& data,
    bool reversed,
    std::vector<double>& out_ts,
    std::vector<Vector3>& out_normals,
    bool record_normals
){
    const size_t n = data.ts.size();
    const bool have_normals = record_normals && data.normals.size() == n;

    out_ts.clear();
    out_normals.clear();

    if (!reversed) {
        out_ts = data.ts;
        if (have_normals) out_normals = data.normals;
        return;
    }

    // Reverse order and flip t-values: t -> 1-t, so the slot stays ascending.
    out_ts.resize(n);
    for (size_t k = 0; k < n; ++k)
        out_ts[k] = 1.0 - data.ts[n - 1 - k];

    if (have_normals) {
        out_normals.resize(n);
        for (size_t k = 0; k < n; ++k)
            out_normals[k] = data.normals[n - 1 - k];
    }
}


// ============================================================================
// SlabEdgeCache
// ============================================================================

void SlabEdgeCache::fill(
    InputQueryHandler& handler,
    const std::vector<std::pair<size_t,size_t>>& edges,
    Map& into, bool use_robust, bool record_normals
){
    into.clear();
    if (edges.empty()) return;

    const size_t n = grid.resolution();

    // The edge list is known up front, so every thread writes its own slots of a
    // pre-sized vector: no insert, no erase, no lock. The map is built afterwards
    // from the results, single-threaded and cheap next to the queries themselves.
    std::vector<EdgeIsect> results(edges.size());
    parallel_for(0, edges.size(), num_threads, [&](size_t k) {
        const auto [i, j] = edges[k];
        const auto ca = grid_node_index_to_coords(n, i);
        const auto cb = grid_node_index_to_coords(n, j);
        const Vector3 pi = grid_node_position(n, true, ca[0], ca[1], ca[2]);
        const Vector3 pj = grid_node_position(n, true, cb[0], cb[1], cb[2]);
        handler.query_edge(i, j, pi, pj, results[k].ts, results[k].normals,
                           use_robust, record_normals);
    });
    num_edge_queries += edges.size();

    // Only edges with crossings are stored; within a loaded plane a miss means
    // "queried, none".
    for (size_t k = 0; k < edges.size(); ++k)
        if (!results[k].ts.empty())
            into.emplace(edge_key(edges[k].first, edges[k].second), std::move(results[k]));
}


void SlabEdgeCache::begin_chunk(
    InputQueryHandler& handler, size_t chunk,
    bool use_robust, bool record_normals
){
    if (loaded && chunk == cur_z) return;

    if (loaded && chunk == cur_z + 1) {
        // Slide: the old upper plane is the new lower one, already queried.
        lo_plane = std::move(hi_plane);
    } else {
        // First chunk, or a jump: rebuild the lower plane from scratch.
        grid.plane_edges(chunk, scratch);
        fill(handler, scratch, lo_plane, use_robust, record_normals);
    }

    grid.cross_edges(chunk, scratch);
    fill(handler, scratch, cross, use_robust, record_normals);

    grid.plane_edges(chunk + 1, scratch);
    fill(handler, scratch, hi_plane, use_robust, record_normals);

    cur_z = chunk;
    loaded = true;
}


void SlabEdgeCache::query_tet(
    InputQueryHandler& /*handler*/,
    const std::array<size_t,4>& tet_indices,
    const std::array<Vector3,4>& /*tet_positions*/,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool /*use_robust*/, bool record_normals
){
    if (!loaded)
        throw std::logic_error("SlabEdgeCache::query_tet before begin_chunk");

    const size_t n = grid.resolution();
    const int zlo = (int)cur_z, zhi = (int)cur_z + 1;

    for (int e = 0; e < 6; ++e) {
        const size_t gi = tet_indices[ALL_TET_PAIRS[e].first];
        const size_t gj = tet_indices[ALL_TET_PAIRS[e].second];

        // Dispatching on the endpoints' planes picks the one map that can hold
        // this edge -- and catches a tet handed to the wrong chunk, which would
        // otherwise read as "no intersections" and silently corrupt the output.
        const int za = grid_node_index_to_coords(n, gi)[2];
        const int zb = grid_node_index_to_coords(n, gj)[2];
        if (za < zlo || za > zhi || zb < zlo || zb > zhi)
            throw std::logic_error("SlabEdgeCache: tet edge outside the loaded window");

        const Map& m = (za == zlo && zb == zlo) ? lo_plane
                     : (za == zhi && zb == zhi) ? hi_plane
                                                : cross;

        auto it = m.find(edge_key(gi, gj));
        if (it == m.end()) {
            edge_isect_ts[e].clear();
            edge_isect_normals[e].clear();
            continue;
        }
        emit_edge_isect(it->second, /*reversed=*/gi > gj,
                        edge_isect_ts[e], edge_isect_normals[e], record_normals);
    }
}


// ============================================================================
// Factories
// ============================================================================

std::unique_ptr<EdgeIsectCache> make_edge_cache(
    QueryCache kind, const TetGridRange& range, int num_threads
){
    switch (kind) {
        case QueryCache::NONE:
            return nullptr;
        case QueryCache::SLAB:
            return std::make_unique<SlabEdgeCache>(range, num_threads);
    }
    return nullptr;
}

std::unique_ptr<EdgeIsectCache> make_edge_cache(
    QueryCache kind, const ExplicitTetRange& /*range*/, int /*num_threads*/
){
    if (kind == QueryCache::SLAB)
        log_warn("query cache 'slab' needs the implicit grid's node planes and is "
                 "ignored for explicit tet meshes, whose intersections are already "
                 "precomputed per edge.");
    return nullptr;
}
