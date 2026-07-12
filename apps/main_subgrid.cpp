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
// negligible effect on the output / non-even-tet count (welding the input is
// what actually matters). A developer can flip this to true if they need it.
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

    std::string out_dir = "./meshes/results/", out_name = "output";
    bool save_output = false;
    if (outputMeshFilename){
        ensure_path_exists(args::get(outputMeshFilename));
        parse_output_filename_and_dir(args::get(outputMeshFilename), out_dir, out_name);
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
        std::cerr << "Error: Provide exactly one of --input (mesh), --inputSDF, or --npz." << std::endl;
        return 1;
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
            query_handler = std::make_unique<SDFQueryHandler>(args::get(inputSDFName), 1e-2f);
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
    std::string filename = "";
    if (save_output){
        std::string out_dir_sub = out_dir + (mod2Flag ? "/mod2" : "/subgrid");
        ensure_path_exists(out_dir_sub + "/dummy.dummy");
        filename = out_dir_sub + "/" + out_name + ".obj";
        std::cout << " Will save output mesh to " << filename << "\n";
    }
    if (TriangleSoup::COMB_MERGE){
        auto [mesh, geo] = makeSurfaceMeshAndGeometry(global_soup.faces, global_soup.vertices);
        mesh->greedilyOrientFaces();
#ifdef HAVE_POLYSCOPE
        if (!noVisFlag) polyscope::registerSurfaceMesh("Subgrid Output Mesh [COMB MERGRED]", geo->inputVertexPositions, mesh->getFaceVertexList())->setSurfaceColor({0.3, 0.6, 0.2})->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
#endif
        if (save_output) writeSurfaceMesh(*mesh, *geo, filename);
    }
    else if (numerical_merge) {
        auto [mesh, geo] = mergeIdenticalVertices(EPS, global_soup.faces, global_soup.vertices);
#ifdef HAVE_POLYSCOPE
        if (!noVisFlag) polyscope::registerSurfaceMesh("Subgrid Output Mesh [NUMERICAL MERGE]", geo->inputVertexPositions, mesh->getFaceVertexList())->setSurfaceColor({0.3, 0.6, 0.2})->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
#endif
        if (save_output) writeSurfaceMesh(*mesh, *geo, filename);
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
        std::cout << ANSI_FG_YELLOW << " [warning] even-sum condition not satisfied; output might not be orientable and might have open holes."
                  << " This is likely due to a non-watertight input or a grid edge-intersection issue. " << "\n"
                  << " \t IF input is watertight, changing the resolution or the perturbation could help." << ANSI_RESET << "\n";
    std::cout << " output: " << global_soup.vertices.size() << " verts, " << global_soup.faces.size() << " faces\n";
    std::cout << " --------------------------------- \n";
    std::cout << " query: "        << result.isect_time << "s"
              << "  construction: " << result.construction_time    << "s"
              << "  overall: "      << overall_time     << "s\n";

#ifdef HAVE_POLYSCOPE
    if (!noVisFlag) polyscope::show();
#endif
    return EXIT_SUCCESS;
}
