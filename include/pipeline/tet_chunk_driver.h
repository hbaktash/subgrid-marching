#pragma once

#include "common/parallel_for.h"
#include "common/triangle_soup.h"
#include "common/utils.h"
#include "query/edge_isect_cache.h"
#include "query/input_query_handler.h"

#include <algorithm>
#include <array>
#include <vector>

// ============================================================================
// Chunked parallel tet loop
// ============================================================================
//
// Steps 1-3 of the pipeline (query, boundary curves, local construction) are
// independent per tet and run across threads. Step 4, folding each local soup
// into the global one, cannot: it threads every vertex signature through one
// shared map. So the driver alternates -- construct a chunk in parallel, merge
// that chunk in flat tet order on the calling thread, move on.
//
// Merging in flat order is what makes the result identical to the serial
// pipeline: per-tet construction depends only on the tet's own positions and
// global indices, never on iteration history, so the local soups are
// order-independent and replaying them in order reproduces the same global mesh.
//
// Chunking is not an optimization, it is a requirement: every local TriangleSoup
// carries its own signature map, so buffering all of them before merging runs
// past a gigabyte at n=128.

// One tet's query results, reused across chunks to avoid churning allocations.
struct TetQuery {
    std::array<std::vector<double>, 6>  ts;
    std::array<std::vector<Vector3>, 6> normals;

    void clear() {
        for (int e = 0; e < 6; ++e) { ts[e].clear(); normals[e].clear(); }
    }
};

// What a worker produces for one tet, consumed by the sequential merge.
struct LocalTetOutput {
    TriangleSoup soup;
    bool active     = false;  // the tet had at least one edge intersection
    bool emitted    = false;  // it produced a soup worth merging
    bool non_even   = false;  // open boundary curves
    bool non_normal = false;  // open or scoop curves
};

// build(flat_index, tet_data, query) -> LocalTetOutput   runs on worker threads
// merge(LocalTetOutput&&)                                runs in order, caller thread
template <typename TetRange, typename BuildFn, typename MergeFn>
void run_tet_chunks(
    const TetRange& range,
    InputQueryHandler& handler,
    EdgeIsectCache* cache,
    int num_threads,
    size_t chunk_size,
    bool use_robust,
    bool record_normals,
    bool show_progress,
    BuildFn&& build,
    MergeFn&& merge
){
    const size_t total = range.total_tet_count();
    if (total == 0) return;

    // A range with no chunk structure (an explicit tet mesh) still needs a bound
    // on how many local soups are alive at once.
    const size_t chunk = chunk_size ? chunk_size : 8192;

    std::vector<TetQuery>      queries;
    std::vector<LocalTetOutput> locals;

    for (size_t lo = 0; lo < total; lo += chunk) {
        const size_t hi = std::min(lo + chunk, total);
        const size_t count = hi - lo;

        if (cache)
            cache->begin_chunk(handler, chunk_size ? lo / chunk_size : 0,
                               use_robust, record_normals);

        if (queries.size() < count) queries.resize(count);
        locals.assign(count, LocalTetOutput{});

        parallel_for(0, count, num_threads, [&](size_t k) {
            const auto tet = range.tet_at(lo + k);
            TetQuery& q = queries[k];
            // Cleared explicitly: not every handler clears its output slots (the
            // SDF path appends), and these buffers are reused across chunks.
            q.clear();
            // The cache's window is fully loaded by begin_chunk, so this only
            // reads it -- no locking needed on either path.
            if (cache)
                cache->query_tet(handler, tet.indices, tet.positions, q.ts, q.normals,
                                 use_robust, record_normals);
            else
                handler.query_intersections(tet.indices, tet.positions, q.ts, q.normals,
                                            use_robust, record_normals);
            locals[k] = build(lo + k, tet, q);
        });

        for (size_t k = 0; k < count; ++k)
            merge(std::move(locals[k]));

        if (show_progress) print_progress(double(hi) / double(total));
    }
}
