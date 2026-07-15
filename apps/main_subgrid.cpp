#include "args.hxx"

#include "common/visual_utils.h"
#include "common/file_IO.h"

#include "assembly/mesh_processing.h"
#include "query/input_query_handler.h"

#include "subgrid_pipeline.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include <chrono>

using namespace geometrycentral;
using namespace geometrycentral::surface;

auto now = std::chrono::high_resolution_clock::now;
template<typename T> using duration = std::chrono::duration<T>;

// Robust ray-intersection queries (FCPW watertight test + iterative all-hits).
// Not exposed on the CLI on purpose: benchmarking showed it ~2x slower with
// negligible effect on the intersection query robustness.
static constexpr bool USE_ROBUST_QUERIES = false;


// ============================================================================
// main
// ============================================================================

void myCallback() {}

int main(int argc, char** argv) {
    auto very_start_time = now();

    // ---- argument parsing ----
    args::ArgumentParser parser("Extract an isosurface from a scalar function defined on a triangle mesh or SDF.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> inputMeshFilename(parser, "inputMesh", "Input mesh file (OBJ, OFF, PLY)", {'i', "input"});
    args::ValueFlag<std::string> inputSDFName(parser, "inputSDF", "Named built-in SDF (e.g. Sphere, Torus); see deps/sdf-dataset", {'s', "inputSDF"});
    args::ValueFlag<std::string> outputMeshFilename(parser, "outputMesh", "Output mesh file (OBJ, OFF, PLY)", {'o', "output"});
    args::ValueFlag<size_t> tetGridResolution(parser, "tetGridResolution", "Grid resolution (e.g., 32, 64, 128 for n^3 grid)", {'r', "tetGridResolution"}, 64);
    args::ValueFlag<double> mergeEPS(parser, "mergeEPS", "Positional merge epsilon; passing this switches from the default combinatorial merge to numerical merge", {"mergeEPS"}, -1);
    args::Flag noMergeFlag(parser, "noMerge", "Disable all vertex merging (output raw triangle soup)", {"noMerge"});
    args::Flag mod2Flag(parser, "mod2", "Use mod2 intersection reduction before subgrid construction", {"mod2"});
    args::ValueFlag<std::string> inputSaveDir(parser, "inputSaveDir", "Save scaled/preprocessed input mesh to this path and exit (mesh mode only)", {"inputSaveDir"});
    args::Flag greedyFlag(parser, "greedy", "Use greedy construction (fan triangulation from boundary curves)", {"greedy"});
    args::ValueFlag<std::string> inputNpzFilename(parser, "inputNpz", "Explicit tet mesh + precomputed intersections (.npz)", {"npz"});
    args::Flag noVisFlag(parser, "noVis", "Disable visualization", {"noViz"});
    args::Flag noProgBar(parser, "noProgBar", "Disable progress bar", {"noPBar"});
    args::Flag listSDFsFlag(parser, "listSDFs", "List the available built-in SDF names and exit", {"listSDFs"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    try {
        if (listSDFsFlag) {
            std::cout << "Available built-in SDFs:\n";
            for (const auto& n : SDFQueryHandler::available_sdf_names())
                std::cout << "  " << n << "\n";
            return EXIT_SUCCESS;
        }

        std::string out_path;
        bool save_output = false;
        if (outputMeshFilename){
            out_path = args::get(outputMeshFilename);
            ensure_path_exists(out_path);
            save_output = true;
        }

#ifdef HAVE_POLYSCOPE
        if (!noVisFlag){
            polyscope::state::userCallback = myCallback;
            initialize_polyscope();
        }
#endif

        // ---- input mode ----
        bool use_sdf = inputSDFName;
        bool use_mesh = inputMeshFilename;
        bool use_npz = inputNpzFilename;
        int input_count = (int)use_sdf + (int)use_mesh + (int)use_npz;
        if (input_count != 1) {
            log_error("provide exactly one of --input (mesh), --inputSDF, or --npz.");
            return EXIT_FAILURE;
        }
        if (use_sdf && !SDFQueryHandler::is_valid_sdf_name(args::get(inputSDFName))) {
            log_error("unknown SDF '" + args::get(inputSDFName) + "'. Use --listSDFs to see available options.");
            return EXIT_FAILURE;
        }
        if (use_mesh) {
            const std::string mesh_path = args::get(inputMeshFilename);
            if (mesh_path.size() >= 4 && mesh_path.compare(mesh_path.size() - 4, 4, ".npz") == 0) {
                log_error("'" + mesh_path + "' is a .npz; pass precomputed intersections with --npz, not -i (which expects a mesh: OBJ/PLY/OFF).");
                return EXIT_FAILURE;
            }
        }

        // ---- merge policy ----
        // Combinatorial merge is the default. --mergeEPS <eps> switches to positional
        // (numerical) merge; --noMerge disables merging entirely (raw triangle soup).
        double EPS = args::get(mergeEPS);
        bool no_merge = noMergeFlag;
        bool numerical_merge = !no_merge && (EPS >= 0);
        bool comb_merge = !no_merge && !numerical_merge;
        TriangleSoup::COMB_MERGE = comb_merge;

        // ---- run pipeline ----
        SubgridPipelineOpts opts;
        opts.mod2 = mod2Flag;
        opts.greedy = greedyFlag;
        opts.use_robust = USE_ROBUST_QUERIES;
        opts.show_progress = !noProgBar;

        SubgridPipelineResult result;

        if (use_npz) {
            result = run_subgrid_pipeline_npz(args::get(inputNpzFilename), opts);
        } else {
            std::unique_ptr<InputQueryHandler> query_handler;
            if (use_sdf) {
                query_handler = std::make_unique<SDFQueryHandler>(args::get(inputSDFName), 1e-3f);
            } else {
                SimplePolygonMesh simpleMesh;
                simpleMesh.readMeshFromFile(args::get(inputMeshFilename), "");
                auto preprocessed = preprocess_input_mesh(simpleMesh);
                if (inputSaveDir) {
                    save_polygon_soup_as_obj(args::get(inputSaveDir), preprocessed.positions, preprocessed.polygons);
                    return EXIT_SUCCESS;
                }
#ifdef HAVE_POLYSCOPE
                if (!noVisFlag)
                    polyscope::registerSurfaceMesh("Input Mesh", preprocessed.positions, preprocessed.polygons)->setEnabled(false);
#endif
                query_handler = std::make_unique<MeshQueryHandler>(preprocessed.positions, preprocessed.polygons);
            }
            result = run_subgrid_pipeline(*query_handler, args::get(tetGridResolution), opts);
        }

        TriangleSoup& global_soup = result.soup;

        // ---- postprocess and save ----
        const std::string& filename = out_path;
        if (TriangleSoup::COMB_MERGE){
            auto [mesh, geo] = makeSurfaceMeshAndGeometry(global_soup.faces, global_soup.vertices);
            mesh->greedilyOrientFaces();
#ifdef HAVE_POLYSCOPE
            if (!noVisFlag) polyscope::registerSurfaceMesh("Subgrid Output Mesh [COMB MERGRED]", geo->inputVertexPositions, mesh->getFaceVertexList())->setSurfaceColor({0.3, 0.6, 0.2})->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
#endif
            if (save_output) writeSurfaceMesh(*mesh, *geo, filename, "obj");
        }
        else if (numerical_merge) {
            auto [mesh, geo] = mergeIdenticalVertices(EPS, global_soup.faces, global_soup.vertices);
#ifdef HAVE_POLYSCOPE
            if (!noVisFlag) polyscope::registerSurfaceMesh("Subgrid Output Mesh [NUMERICAL MERGE]", geo->inputVertexPositions, mesh->getFaceVertexList())->setSurfaceColor({0.3, 0.6, 0.2})->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
#endif
            if (save_output) writeSurfaceMesh(*mesh, *geo, filename, "obj");
        }
        else {
#ifdef HAVE_POLYSCOPE
            if (!noVisFlag) polyscope::registerSurfaceMesh("Subgrid Output Mesh [NO MERGE]", global_soup.vertices, global_soup.faces)->setSurfaceColor({0.3, 0.6, 0.2})->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
#endif
            if (save_output) save_polygon_soup_as_obj(filename, global_soup.vertices, global_soup.faces);
        }

        double overall_time = duration<double>(now() - very_start_time).count();
        std::cout.precision(4);
        std::cout << " non-zero tets:   " << result.non_zero_tets << "/" << result.total_tets << "\n";
        std::cout << " non-normal tets: " << result.non_normal_tets << "/" << result.non_zero_tets << "\n";
        std::cout << " non-even tets:   " << result.non_even_tets  << "/" << result.non_zero_tets << "\n";
        if (result.non_even_tets > 0)
            log_warn("even-sum condition not satisfied; output might not be orientable and might have open holes."
                     " Likely a non-watertight input or a grid edge-intersection issue."
                     " If the input is watertight, changing the resolution or the perturbation could help.");
        std::cout << " output: " << global_soup.vertices.size() << " verts, " << global_soup.faces.size() << " faces\n";
        if (save_output) std::cout << " saved output mesh to " << filename << "\n";
        std::cout << " --------------------------------- \n";
        std::cout << " query: "        << result.isect_time << "s"
                  << "  construction: " << result.construction_time    << "s"
                  << "  overall: "      << overall_time     << "s\n";

#ifdef HAVE_POLYSCOPE
        if (!noVisFlag) polyscope::show();
#endif
    } catch (const std::exception& e) {
        log_error(e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
