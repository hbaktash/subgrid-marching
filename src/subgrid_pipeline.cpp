#include "subgrid_pipeline.h"
#include "grid/grid_iterators.h"
#include "grid/explicit_tet_range.h"
#include "query/precomputed_query_handler.h"
#include "query/intersection_query.h"
#include "query/mod2_reduction.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/primal_reconstruction.h"
#include "common/npz_reader.h"
#include "common/utils.h"

#include "query/edge_isect_cache.h"
#include "pipeline/tet_chunk_driver.h"

#include <chrono>
#include <limits>
#include <numeric>

namespace {
auto clock_now = std::chrono::high_resolution_clock::now;
template<typename T> using duration = std::chrono::duration<T>;
}


// Everything a worker does for one tet: steps 2 and 3 of the pipeline. Pure in
// its inputs (the tet's own positions, global indices, and edge intersections),
// which is exactly why the chunks can run out of order and still merge to the
// same mesh.
static LocalTetOutput build_primal_tet(
    const std::array<Vector3,4>& tet_positions,
    const std::array<size_t,4>& tet_indices,
    std::array<std::vector<double>,6>& edge_isect_ts,
    const SubgridPipelineOpts& opts
){
    LocalTetOutput out;

    auto edge_intersection_counts = counts_from_isect_ts(edge_isect_ts);
    if (std::accumulate(edge_intersection_counts.begin(), edge_intersection_counts.end(), 0) == 0)
        return out;
    out.active = true;

    if (opts.mod2 && apply_mod2_reduction(edge_isect_ts, edge_intersection_counts))
        return out;

    EdgeOccupations non_normal_edge_occupations;
    auto [open_curves, scoop_curves, normal_curves] = boundary_comb_curves(
        tet_indices, edge_intersection_counts, non_normal_edge_occupations, opts.greedy
    );
    out.non_even   = !open_curves.empty();
    out.non_normal = !(open_curves.empty() && scoop_curves.empty());

    if (opts.greedy) {
        out.soup = subgrid_greedy(
            tet_positions, tet_indices, edge_isect_ts,
            scoop_curves, normal_curves,
            opts.scoop_mid_vertices, opts.scoop_bulge
        );
    } else {
        out.soup = subgrid_surface(
            tet_positions, tet_indices, edge_isect_ts,
            open_curves, scoop_curves, non_normal_edge_occupations,
            opts.scoop_mid_vertices, opts.scoop_bulge, opts.verbose
        );
    }
    out.emitted = true;
    return out;
}


// Threaded variant of the loop below. Same steps in the same order; only the
// scheduling differs, and the merge still runs sequentially in tet order.
template<typename TetRange>
static SubgridPipelineResult run_subgrid_impl_mt(
    TetRange& tet_range, InputQueryHandler& handler, const SubgridPipelineOpts& opts
){
    SubgridPipelineResult result;
    result.total_tets = tet_range.total_tet_count();

    auto cache = make_edge_cache(opts.query_cache, tet_range, opts.num_threads);
    const size_t chunk_size = cache ? cache_chunk_size(tet_range) : 0;

    // The handler's diagnostic counters are plain shared state; nothing writes
    // them while threads are running.
    // Diagnostic counters are plain shared state, so they stay off while threads run.
    QueryModeGuard query_mode(handler, opts.canonical_queries, /*stats=*/false);

    const auto wall_start = clock_now();
    double merge_time = 0.0;

    run_tet_chunks(
        tet_range, handler, cache.get(), opts.num_threads, chunk_size,
        opts.use_robust, /*record_normals=*/false, opts.show_progress,
        [&](size_t, const auto& tet, TetQuery& q) {
            return build_primal_tet(tet.positions, tet.indices, q.ts, opts);
        },
        [&](LocalTetOutput&& local) {
            const auto t0 = clock_now();
            if (local.active)     result.non_zero_tets++;
            if (local.non_even)   result.non_even_tets++;
            if (local.non_normal) result.non_normal_tets++;
            if (local.emitted)    result.soup.add_local_soup(local.soup);
            merge_time += duration<double>(clock_now() - t0).count();
        }
    );

    if (opts.show_progress) std::cout << std::endl;

    // Wall-clock, not summed per-thread time: with threads running the two
    // phases overlap, so a per-phase breakdown would add up to more than elapsed.
    result.construction_time = merge_time;
    result.isect_time = duration<double>(clock_now() - wall_start).count() - merge_time;
    return result;
}


template<typename TetRange>
static SubgridPipelineResult run_subgrid_impl(TetRange& tet_range, InputQueryHandler& handler, const SubgridPipelineOpts& opts) {
    if (opts.num_threads != 1)
        return run_subgrid_impl_mt(tet_range, handler, opts);

    SubgridPipelineResult result;
    result.total_tets = tet_range.total_tet_count();

    // Optional edge intersection reuse. Chunks (cube layers, for the grid) are
    // visited in order, which is what lets the slab window slide.
    auto cache = make_edge_cache(opts.query_cache, tet_range);
    const size_t chunk_size = cache ? cache_chunk_size(tet_range) : 0;
    size_t cur_chunk = std::numeric_limits<size_t>::max();
    QueryModeGuard query_mode(handler, opts.canonical_queries, handler.collect_query_stats);

    auto start_time = clock_now(), end_time = clock_now();
    size_t processed_tets = 0;
    for (const auto& tet_data : tet_range) {
        if (opts.show_progress && (processed_tets % 10000 == 0))
            print_progress(double(processed_tets) / double(result.total_tets));
        const size_t tet_ordinal = processed_tets++;

        std::array<Vector3,4> tet_positions = tet_data.positions;
        std::array<size_t,4>  tet_indices   = tet_data.indices;

        // intersection query
        start_time = clock_now();
        std::array<std::vector<double>,6>  edge_isect_ts;
        std::array<std::vector<Vector3>,6> edge_isect_normals;
        // primal pipeline never uses intersection normals, so skip recording them
        if (cache) {
            const size_t chunk = chunk_size ? tet_ordinal / chunk_size : 0;
            if (chunk != cur_chunk) {
                cache->begin_chunk(handler, chunk, opts.use_robust, /*record_normals=*/false);
                cur_chunk = chunk;
            }
            cache->query_tet(handler, tet_indices, tet_positions, edge_isect_ts, edge_isect_normals,
                             opts.use_robust, /*record_normals=*/false);
        } else {
            handler.query_intersections(tet_indices, tet_positions, edge_isect_ts, edge_isect_normals,
                                        opts.use_robust, /*recordNormals=*/false);
        }
        end_time = clock_now();
        result.isect_time += duration<double>(end_time - start_time).count();

        auto edge_intersection_counts = counts_from_isect_ts(edge_isect_ts);
        bool empty_tet = (std::accumulate(edge_intersection_counts.begin(), edge_intersection_counts.end(), 0) == 0);
        // The per-tet query-count breakdown has no meaning once edges are served
        // from the cache, so it is only collected on the direct path.
        if (handler.is_sdf() && !cache)
            handler.update_global_query_count_map(tet_indices, !empty_tet);
        if (empty_tet) continue;
        result.non_zero_tets++;

        // subgrid construction
        start_time = clock_now();

        if (opts.mod2 && apply_mod2_reduction(edge_isect_ts, edge_intersection_counts))
            continue;

        EdgeOccupations non_normal_edge_occupations;
        auto [open_curves, scoop_curves, normal_curves] = boundary_comb_curves(
            tet_indices, edge_intersection_counts, non_normal_edge_occupations, opts.greedy
        );
        result.non_even_tets   += open_curves.empty() ? 0 : 1;
        result.non_normal_tets += (open_curves.empty() && scoop_curves.empty()) ? 0 : 1;

        TriangleSoup local_soup;
        if (opts.greedy) {
            local_soup = subgrid_greedy(
                tet_positions, tet_indices, edge_isect_ts,
                scoop_curves, normal_curves,
                opts.scoop_mid_vertices, opts.scoop_bulge
            );
        } else {
            local_soup = subgrid_surface(
                tet_positions, tet_indices, edge_isect_ts,
                open_curves, scoop_curves, non_normal_edge_occupations,
                opts.scoop_mid_vertices, opts.scoop_bulge, opts.verbose
            );
        }
        result.soup.add_local_soup(local_soup);

        end_time = clock_now();
        result.construction_time += duration<double>(end_time - start_time).count();
    }
    if (opts.show_progress) { print_progress(1.0); std::cout << std::endl; }
    return result;
}

SubgridPipelineResult run_subgrid_pipeline(InputQueryHandler& handler, size_t resolution, const SubgridPipelineOpts& opts) {
    TetGridRange tet_range(resolution, true);
    return run_subgrid_impl(tet_range, handler, opts);
}

SubgridPipelineResult run_subgrid_pipeline_npz(const std::string& npz_path, const SubgridPipelineOpts& opts) {
    auto npz = load_npz(npz_path);
    ExplicitTetRange tet_range(npz["vertices"], npz["tets"]);
    const NpyArray* normals_ptr = npz.has("isect_normals") ? &npz["isect_normals"] : nullptr;
    PrecomputedQueryHandler handler(npz["edges"], npz["isect_offsets"], npz["isect_ts"], normals_ptr);
    return run_subgrid_impl(tet_range, handler, opts);
}
