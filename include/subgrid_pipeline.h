#pragma once

#include "subgrid_MT/cv_interpolation.h"
#include "query/input_query_handler.h"
#include <cstddef>
#include <string>

// Options for the primal subgrid pipeline. Defaults match the canonical single-tet
// construction (see docs/construction_policy.md).
struct SubgridPipelineOpts {
    bool mod2 = false;                 // reduce each edge count mod 2 before construction
    bool greedy = false;               // greedy fan construction instead of subgrid_surface
    bool use_robust = false;           // robust ray-intersection queries (mesh input)
    bool show_progress = false;        // print a progress bar during the tet loop
    bool verbose = false;
    double scoop_bulge = 0.001;       // Steiner mid-scoop-vertex bulge magnitude
    bool scoop_mid_vertices = true;    // simplicial-embedding Steiner insertion
};

struct SubgridPipelineResult {
    TriangleSoup soup;
    int non_zero_tets = 0;
    int non_normal_tets = 0;
    int non_even_tets = 0;  // tets that had open curves (odd face-sum somewhere)
    size_t total_tets = 0;
    double isect_time = 0.0;
    double construction_time = 0.0;
};

// Run the primal subgrid pipeline over a full tet grid.
// non_even_tets == 0 iff every active tet is even-sum, in which case the merged
// output is guaranteed to be orientable.
// Combinatorial merging follows the global TriangleSoup::COMB_MERGE flag; set it
// before calling if exact merging is desired.
SubgridPipelineResult run_subgrid_pipeline(InputQueryHandler& handler, size_t resolution, const SubgridPipelineOpts& opts = {});

// Run the primal subgrid pipeline over an explicit tet mesh with precomputed intersections.
SubgridPipelineResult run_subgrid_pipeline_npz(const std::string& npz_path, const SubgridPipelineOpts& opts = {});
