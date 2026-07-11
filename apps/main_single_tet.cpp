#include "args.hxx"

#include "common/visual_utils.h"
#include "grid/mesh_factory.h"


// new stuff
#include "subgrid_MT/primal_reconstruction.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/cv_interpolation.h"
#include "subgrid_MT/evensum_construction.h"

#include "common/file_IO.h"

#include <filesystem>
#include <fstream>
#include <map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

int ne_01 = 1, ne_02 = 2, ne_03 = 3,
    ne_12 = 7, ne_13 = 2, 
    ne_23 = 1;
int ne_x1 = 0, ne_x2 = 0, ne_x3 = 0;

bool scoop_mid_vertices = true;
float scoop_mid_vertex_bulge = 0.01;

int corner_cuts[4] = {0, 0, 0, 0};
int diagonal_cuts[3] = {0, 4, 0};

int max_nc = 10;

bool allow_intersections = false,
     only_recursive_build = false,
     always_show_other_tet = false;

// When loaded from a config file, overrides uniform t-values.
// Cleared on any slider change, reverting to uniform.
std::array<std::vector<double>,6> loaded_isect_ts;
bool use_loaded_ts = false;

/// output soup
std::vector<std::vector<size_t>> output_face_indices;
std::vector<Vector3> output_vertex_positions;
std::string out_path = "";
std::string stats_out_path = ""; // empty = stdout only; unset when stats_enabled=false
bool stats_enabled = false;
bool noViz = false;

using std::vector;
using std::array;

std::array<std::vector<Vector3>,6> 
generate_dummy_normals(
    const std::array<std::vector<double>,6> &edge_isect_ts
){
    std::array<std::vector<Vector3>,6> edge_isect_normals;
    Vector3 default_normal{0., 0., 0.};
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        int n_ij = edge_isect_ts[eIdx].size();
        edge_isect_normals[eIdx] = std::vector<Vector3>(n_ij, default_normal);
    }
    return edge_isect_normals;
}

array<vector<double>,6> 
generate_uniform_isect_ts(
    EdgeInts edge_ints
){
    array<vector<double>,6> edge_isect_ts;
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        int i = ALL_TET_PAIRS[eIdx].first,
            j = ALL_TET_PAIRS[eIdx].second;
        int nij = edge_ints[{i,j}];
        for (int k = 0; k < nij; ++k) {
            edge_isect_ts[eIdx].push_back(double(k + 1) / double(nij + 1));
        }
    }
    return edge_isect_ts;
}

void remove_other_tet_structures(){
    polyscope::removeStructure("Other Tet");
    polyscope::removeStructure("Subgrid Surface_xTet");
    polyscope::removeStructure("Open Curves_xTet");
    polyscope::removeStructure("Scoop Curves_xTet");
    polyscope::removeStructure("Normal Curves_xTet");
    for (size_t e = 0; e < 6; ++e){
        std::string name = "Edge " + std::to_string(ALL_TET_PAIRS[e].first) + "-" + std::to_string(ALL_TET_PAIRS[e].second) + " Isects_xTet";
        if (polyscope::hasPointCloud(name))
            polyscope::removePointCloud(name);
    }
}

void 
embed_and_visualize_boundary_curves(
    const vector<CombFace> &open_comb_curves,
    const vector<CombFace> &scoop_comb_curves,
    const vector<CombFace> &normal_comb_curves,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts,
    double bulge,
    std::string name_suffix = ""
){
    vector<vector<Vector3>> open_curve_positions, scoop_curve_positions, normal_curve_positions;
    for (const auto &comb_curve: open_comb_curves) {
        open_curve_positions.push_back(comb_curve_to_vec3_curve(comb_curve, tet_positions, edge_isect_ts, bulge));
    }
    for (const auto &comb_curve: scoop_comb_curves) {
        scoop_curve_positions.push_back(comb_curve_to_vec3_curve(comb_curve, tet_positions, edge_isect_ts, bulge));
    }
    for (const auto &comb_curve: normal_comb_curves) {
        normal_curve_positions.push_back(comb_curve_to_vec3_curve(comb_curve, tet_positions, edge_isect_ts, bulge));
    }
    visualizeBoundaryCurves(open_curve_positions, false, glm::vec3(1.0, 0.2, 0.2), "Open Curves" + name_suffix, 0.004);
    visualizeBoundaryCurves(scoop_curve_positions, true, glm::vec3(0.4, 1.0, 0.4), "Scoop Curves" + name_suffix, 0.004);
    visualizeBoundaryCurves(normal_curve_positions, true, glm::vec3(0.2, 0.2, 1.0), "Normal Curves" + name_suffix, 0.004);
    // polyscope::show();
}

EdgeInts 
get_edge_ints_from_user_input(){
    EdgeInts edge_ints;
    edge_ints.set(0,1, ne_01);
    edge_ints.set(0,2, ne_02);
    edge_ints.set(0,3, ne_03);
    edge_ints.set(1,2, ne_12);
    edge_ints.set(1,3, ne_13);
    edge_ints.set(2,3, ne_23);
    return edge_ints;
}

EdgeInts 
get_other_tet_edge_ints_from_user_input(){
    EdgeInts edge_ints;
    edge_ints.set(0,1, ne_x1);
    edge_ints.set(0,2, ne_x2);
    edge_ints.set(0,3, ne_x3);
    edge_ints.set(1,2, ne_12);
    edge_ints.set(1,3, ne_13);
    edge_ints.set(2,3, ne_23);
    return edge_ints;
}

void
print_curve_stats(
    const vector<CombFace> &open_curves,
    const vector<CombFace> &scoop_curves,
    const vector<CombFace> &normal_curves,
    std::ostream &out
){
    auto print_size_hist = [&](const vector<CombFace> &curves, const std::string &label){
        out << label << ":\n";
        if (curves.empty()){ out << "  (none)\n"; return; }
        std::map<int,int> hist;
        for (const auto &c: curves) hist[(int)c.vertices.size()]++;
        for (const auto &[sz, cnt]: hist) out << "  size " << sz << ": " << cnt << "\n";
    };

    print_size_hist(open_curves,   "open curves");
    print_size_hist(scoop_curves,  "scoop curves");
    print_size_hist(normal_curves, "normal curves");

    out << "scoop curve types:\n";
    int contractible = 0, corner_type = 0, diag_type = 0;
    for (const auto &sc: scoop_curves){
        int sep = (int)separated_vertices(combface_to_edgeints(sc)).size();
        if      (sep == 0 || sep == 4) contractible++;
        else if (sep == 1 || sep == 3) corner_type++;
        else if (sep == 2)             diag_type++;
    }
    out << "  contractible: " << contractible << "\n";
    out << "  corner-type:  " << corner_type  << "\n";
    out << "  diag-type:    " << diag_type     << "\n";
}

void
log_curve_stats(
    const vector<CombFace> &open_curves,
    const vector<CombFace> &scoop_curves,
    const vector<CombFace> &normal_curves
){
    if (!stats_enabled) return;

    print_curve_stats(open_curves, scoop_curves, normal_curves, std::cout);
    if (stats_out_path.empty()) return;

    std::string stem, parent, extension;
    get_stem_parent_extension(stats_out_path, stem, parent, extension);
    std::string stats_name = "stats_" + stem + "_eints_" +
        std::to_string(ne_01) + "_" + std::to_string(ne_02) + "_" +
        std::to_string(ne_03) + "_" + std::to_string(ne_12) + "_" +
        std::to_string(ne_13) + "_" + std::to_string(ne_23);
    std::string resolved_path = parent + "/" + stats_name + ".txt";
    std::ofstream f(resolved_path);
    if (f) {
        print_curve_stats(open_curves, scoop_curves, normal_curves, f);
        std::cout << "Saved stats to " << resolved_path << "\n";
    } else {
        std::cerr << "Failed to open stats file: " << resolved_path << "\n";
    }
}

void
new_reconstruction_from_single_tet_edge_ints(bool other_tet = false){
    EdgeInts edge_ints = get_edge_ints_from_user_input();
    array<Vector3,4> tet_positions = getTetVertexPositions();
    if (other_tet){
        edge_ints = get_other_tet_edge_ints_from_user_input();
        tet_positions[0] *= -1.;
        if (!noViz)
            polyscope::registerSurfaceMesh("Other Tet", tet_positions, vector<vector<size_t>>{{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}})->setSurfaceColor({0.1, 0.5, 1.})->setBackFacePolicy(polyscope::BackFacePolicy::Identical)->setTransparency(0.5f)->setEnabled(true);
    }
    // Use loaded t-values if available; otherwise generate uniform
    array<vector<double>,6> edge_isect_ts;
    if (use_loaded_ts && !other_tet)
        edge_isect_ts = loaded_isect_ts;
    else
        edge_isect_ts = generate_uniform_isect_ts(edge_ints);

    std::cout << "Edge intersection counts: ";
    for (const auto &ts: edge_isect_ts) {
        std::cout << ts.size() << " ";
    }
    std::cout << "\n";

    if (!noViz) {
        // scaling the tet positions to properly visualize smeared (non-normal) polygons
        for (Vector3 &pos: tet_positions) pos *= 1.001;
        visualizeIntersections(edge_isect_ts, tet_positions, other_tet ? "_xTet" : "");
    }

    // cv curves
    std::cout << "Computing boundary curves for edge ints: ";
    EdgeOccupations non_normal_edge_occupations;
    auto [open_curves, scoop_curves, normal_curves] =
    boundary_comb_curves(
        {0,1,2,3}, edge_ints.to_array(), non_normal_edge_occupations, true
    );

    std::cout << "Showing boundary curves: ";
    if (!noViz)
        embed_and_visualize_boundary_curves(open_curves, scoop_curves, normal_curves,
                tet_positions, edge_isect_ts, (double)scoop_mid_vertex_bulge, other_tet ? "_xTet" : "");

    if (!other_tet) log_curve_stats(open_curves, scoop_curves, normal_curves);

    bool verbose = true;
    std::cout << "Constructing subgrid surface: ";
    auto local_soup = subgrid_surface(tet_positions, {0,1,2,3}, edge_isect_ts, open_curves, scoop_curves, non_normal_edge_occupations,
                        scoop_mid_vertices, (double)scoop_mid_vertex_bulge, verbose);
    auto [face_inds, vertex_positions] = local_soup.get_soup();
    if (!other_tet){
        output_face_indices    = face_inds;
        output_vertex_positions = vertex_positions;
    }
    if (!noViz) {
        std::string title = std::string("Subgrid Surface") + (other_tet ? "_xTet" : "");
        auto* ps_mesh = polyscope::registerSurfaceMesh(title, vertex_positions, face_inds);
        ps_mesh->setSurfaceColor(glm::vec3(0.8, 0.2, 0.2))->setEnabled(true);

    }
}

EdgeInts build_edge_ints(){
    CornerCuts cc;
    cc[0] = corner_cuts[0];
    cc[1] = corner_cuts[1];
    cc[2] = corner_cuts[2];
    cc[3] = corner_cuts[3];
    DiagonalCuts dc;
    dc.set(0,1, diagonal_cuts[0]);
    dc.set(0,2, diagonal_cuts[1]);
    dc.set(0,3, diagonal_cuts[2]);
    EdgeInts edge_ints = corner_and_diagonal_cuts_to_edge_ints(cc, dc);
    ne_01 = edge_ints.get(0,1);
    ne_02 = edge_ints.get(0,2);
    ne_03 = edge_ints.get(0,3);
    ne_12 = edge_ints.get(1,2);
    ne_13 = edge_ints.get(1,3);
    ne_23 = edge_ints.get(2,3);
    return edge_ints;
}

EdgeInts get_first_edge_ints(){
    EdgeInts edge_ints;
    edge_ints.set(0,1, ne_01);
    edge_ints.set(0,2, ne_02);
    edge_ints.set(0,3, ne_03);
    edge_ints.set(1,2, ne_12);
    edge_ints.set(1,3, ne_13);
    edge_ints.set(2,3, ne_23);
    return edge_ints;
}

EdgeInts get_second_edge_ints(){
    EdgeInts edge_ints;
    edge_ints.set(0,1, ne_x1);
    edge_ints.set(0,2, ne_x2);
    edge_ints.set(0,3, ne_x3);
    edge_ints.set(1,2, ne_12);
    edge_ints.set(1,3, ne_13);
    edge_ints.set(2,3, ne_23);
    return edge_ints;
}


void myCallback() {
    if (ImGui::SliderInt4("corner cuts", corner_cuts, 0, max_nc) ||
            ImGui::SliderInt3("diagonal cuts", diagonal_cuts, 0, max_nc)) {
        use_loaded_ts = false;
        EdgeInts edge_ints = build_edge_ints();
        new_reconstruction_from_single_tet_edge_ints();
    }
    ImGui::Checkbox("allow intersections", &allow_intersections);
    ImGui::Checkbox("only recursive build", &only_recursive_build);
    if (ImGui::Checkbox("scoop mid vertices", &scoop_mid_vertices)){
        new_reconstruction_from_single_tet_edge_ints();
        if (always_show_other_tet || ne_x1 > 0 || ne_x2 > 0 || ne_x3 > 0)
            new_reconstruction_from_single_tet_edge_ints(true);
    }
    if (ImGui::Checkbox("comb merge", &TriangleSoup::COMB_MERGE)){
        new_reconstruction_from_single_tet_edge_ints();
        if (always_show_other_tet || ne_x1 > 0 || ne_x2 > 0 || ne_x3 > 0)
            new_reconstruction_from_single_tet_edge_ints(true);
    }
    if (ImGui::SliderFloat("scoop bulge", &scoop_mid_vertex_bulge, 0.0, 0.2)){
        new_reconstruction_from_single_tet_edge_ints();
        if (always_show_other_tet || ne_x1 > 0 || ne_x2 > 0 || ne_x3 > 0)
            new_reconstruction_from_single_tet_edge_ints(true);
    }
    if (ImGui::SliderInt("edge 01", &ne_01, 0, max_nc) ||
        ImGui::SliderInt("edge 02", &ne_02, 0, max_nc) ||
        ImGui::SliderInt("edge 03", &ne_03, 0, max_nc)){
        use_loaded_ts = false;
        new_reconstruction_from_single_tet_edge_ints();
        if (always_show_other_tet)
            new_reconstruction_from_single_tet_edge_ints(true);
    }
    if (ImGui::SliderInt("edge 12", &ne_12, 0, max_nc) ||
        ImGui::SliderInt("edge 13", &ne_13, 0, max_nc) ||
        ImGui::SliderInt("edge 23", &ne_23, 0, max_nc)){
        use_loaded_ts = false;
        new_reconstruction_from_single_tet_edge_ints();
        if (ne_x1 > 0 || ne_x2 > 0 || ne_x3 > 0 || always_show_other_tet)
            new_reconstruction_from_single_tet_edge_ints(true);
    }
    if (ImGui::Checkbox("always show other tet", &always_show_other_tet)){
        if (always_show_other_tet)
            new_reconstruction_from_single_tet_edge_ints(true);
        else if (ne_x1 == 0 && ne_x2 == 0 && ne_x3 == 0)
            remove_other_tet_structures();
    }
    if (ImGui::SliderInt("edge x1", &ne_x1, 0, max_nc) ||
            ImGui::SliderInt("edge x2", &ne_x2, 0, max_nc) ||
            ImGui::SliderInt("edge x3", &ne_x3, 0, max_nc)) {
        new_reconstruction_from_single_tet_edge_ints(true);
        if (ne_x1 == 0 && ne_x2 == 0 && ne_x3 == 0 && !always_show_other_tet){
            remove_other_tet_structures();
        }
    }
    if (ImGui::Button("save soup as OBJ")) {
        if (out_path == "") {
            std::cerr << "No output filename specified; not saving." << std::endl;
        }
        else {
            std::filesystem::path filepath(out_path);
            // get stem, parent and extension
            std::string stem, parent, extension;
            get_stem_parent_extension(out_path, stem, parent, extension);
            // new stem with edge ints
            std::string new_stem = stem + "_eints_" +
                std::to_string(ne_01) + "_" + std::to_string(ne_02) + "_" +
                std::to_string(ne_03) + "_" + std::to_string(ne_12) + "_" +
                std::to_string(ne_13) + "_" + std::to_string(ne_23);
            std::string new_filename = parent + "/" + new_stem + extension;
            save_polygon_soup_as_obj(new_filename, output_vertex_positions, output_face_indices);
            std::cout << "Saved output soup to " << new_filename << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("Interactive normal surface construction for a single tetrahedron.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> outputMeshFilename(parser, "outputMesh", "Output mesh file (OBJ, OFF, PLY)", {'o', "output"});
    args::NargsValueFlag<int> edgeIntsFlag(parser, "edge_ints",
        "Initial edge intersection counts: n01 n02 n03 n12 n13 n23",
        {'e', "edge_ints"}, {6, 6});
    args::ImplicitValueFlag<std::string> statsFlag(parser, "stats",
        "Print curve statistics. Optionally provide dir/name.txt to also write to dir/stats_name_eints_*.txt.",
        {"stats"}, std::string(""), std::string(""));
    args::Flag noVizFlag(parser, "noViz", "Skip Polyscope visualization (headless stats collection)", {"noViz"});
    args::ValueFlag<std::string> loadFileFlag(parser, "config_file", "Load a .txt with edge_ints + t-values", {"load-file"});

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
    if (outputMeshFilename){
        out_path = args::get(outputMeshFilename);
    }
    if (edgeIntsFlag){
        auto vals = args::get(edgeIntsFlag);
        ne_01 = vals[0]; ne_02 = vals[1]; ne_03 = vals[2];
        ne_12 = vals[3]; ne_13 = vals[4]; ne_23 = vals[5];
    }
    if (statsFlag && args::get(statsFlag) != "__absent__"){
        stats_enabled = true;
        stats_out_path = args::get(statsFlag); // empty string = stdout only
    }
    noViz = noVizFlag;

    if (loadFileFlag) {
        TetConfigFileData data;
        if (load_tet_config_file(args::get(loadFileFlag), data)) {
            ne_01 = data.edge_ints[0]; ne_02 = data.edge_ints[1]; ne_03 = data.edge_ints[2];
            ne_12 = data.edge_ints[3]; ne_13 = data.edge_ints[4]; ne_23 = data.edge_ints[5];
            loaded_isect_ts = data.edge_isect_ts;
            use_loaded_ts = true;
            std::cout << "Loaded config file: " << args::get(loadFileFlag) << "\n";
        } else {
            std::cerr << "Failed to load config file: " << args::get(loadFileFlag) << "\n";
        }
    }

    if (!noViz) {
        polyscope::state::userCallback = myCallback;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        polyscope::options::ssaaFactor = 2;
        polyscope::init();
        polyscope::registerSurfaceMesh("Tet Mesh", getTetVertexPositions(), getTetFaceIndices())->setSurfaceColor({0.1, 0.5, 1.})->setTransparency(0.5)->setEnabled(true);
    }

    new_reconstruction_from_single_tet_edge_ints();

    if (!noViz)
        polyscope::show();
    return EXIT_SUCCESS;
}
