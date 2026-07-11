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

#include <chrono>
#include <numeric>

namespace {
auto clock_now = std::chrono::high_resolution_clock::now;
template<typename T> using duration = std::chrono::duration<T>;
}


template<typename TetRange>
static SubgridPipelineResult run_subgrid_impl(TetRange& tet_range, InputQueryHandler& handler, const SubgridPipelineOpts& opts) {
    SubgridPipelineResult result;
    result.total_tets = tet_range.total_tet_count();

    auto start_time = clock_now(), end_time = clock_now();
    for (const auto& tet_data : tet_range) {
        if ((tet_data.indices[0] % 10000 == 0) && opts.show_progress)
            print_progress(double(tet_data.indices[0] * 5) / double(result.total_tets));

        std::array<Vector3,4> tet_positions = tet_data.positions;
        std::array<size_t,4>  tet_indices   = tet_data.indices;

        // intersection query
        start_time = clock_now();
        std::array<std::vector<double>,6>  edge_isect_ts;
        std::array<std::vector<Vector3>,6> edge_isect_normals;
        // primal pipeline never uses intersection normals, so skip recording them
        handler.query_intersections(tet_indices, tet_positions, edge_isect_ts, edge_isect_normals,
                                    opts.use_robust, /*recordNormals=*/false);
        end_time = clock_now();
        result.isect_time += duration<double>(end_time - start_time).count();

        auto edge_intersection_counts = counts_from_isect_ts(edge_isect_ts);
        bool empty_tet = (std::accumulate(edge_intersection_counts.begin(), edge_intersection_counts.end(), 0) == 0);
        if (handler.is_sdf())
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
    if (opts.show_progress) print_progress(1.0);
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
