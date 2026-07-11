#pragma once

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
};

struct DualSubgridPipelineResult {
    TriangleSoup soup;  // vertices, faces, faces_per_edge, dual_positions all inside
    int non_zero_tets = 0;
    int non_normal_tets = 0;
    int non_even_tets = 0;
    size_t total_tets = 0;
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
