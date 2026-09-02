#pragma once

#include "query/edge_isect_cache.h"   // QueryCache
#include "query/input_query_handler.h"
#include "assembly/mesh_processing.h"
#include "subgrid_MT/dual/qef_solver.h"
#include "subgrid_MT/dual/dual_construction.h"
#include "common/triangle_soup.h"
#include <cstddef>
#include <string>

// Options for the dual subgrid pipeline.
struct DualSubgridPipelineOpts {
    bool mod2 = false;             // reduce each edge count mod 2 before construction
    bool use_robust = false;       // robust ray-intersection queries (mesh input)
    bool show_progress = false;    // print a progress bar during the tet loop
    double reg_alpha = 0.1;        // QEF regularization weight
    bool project_duals = false;    // clip dual points back inside the tet
    bool no_normal = false;        // ignore normals: place each dual point at the
                                   // boundary-polygon centroid instead of solving a QEF
    QueryCache query_cache = QueryCache::NONE;  // edge intersection reuse
    // Query every edge in the canonical min->max direction and flip for tets that
    // traverse it the other way, so tets sharing an edge always agree on where --
    // and how many times -- the surface crosses it. Costs no extra queries, and a
    // query cache is always canonical, so this is what makes cached and uncached
    // runs agree bit for bit. See InputQueryHandler::canonical_edge_queries.
    bool canonical_queries = true;
    int num_threads = 1;           // worker threads for the tet loop; 0 = all cores
};

struct DualSubgridPipelineResult {
    TriangleSoup soup;  // vertices, faces, faces_per_edge, dual_positions all inside
    int non_zero_tets = 0;
    int non_normal_tets = 0;
    int non_even_tets = 0;
    size_t total_tets = 0;
    // Serial runs (num_threads == 1) split these as the names say. With threads
    // the query and construction phases overlap, so summing per-thread times
    // would exceed elapsed: isect_time then holds the wall-clock time of the
    // parallel query+construction phases and construction_time the sequential
    // merge. Together they still add up to the tet loop's elapsed time.
    double isect_time = 0.0;
    double construction_time = 0.0;
};

DualSubgridPipelineResult run_dual_subgrid_pipeline(
    InputQueryHandler& handler,
    size_t resolution,
    const DualSubgridPipelineOpts& opts = {}
);

DualSubgridPipelineResult run_dual_subgrid_pipeline_npz(
    const std::string& npz_path,
    const DualSubgridPipelineOpts& opts = {}
);
