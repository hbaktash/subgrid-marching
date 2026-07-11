#include "args.hxx"

#include "common/visual_utils.h"
#include "subgrid_MT/primal_reconstruction.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/cv_interpolation.h"
#include "assembly/mesh_processing.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include <random>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ── ring geometry ────────────────────────────────────────────────────────────
// 4 tets sharing edge {0,1}, ring vertices 2,3,4,5 around it.
// Tets: {0,1,2,3}, {0,1,3,4}, {0,1,4,5}, {0,1,5,2}

static const std::array<Vector3,6> RING_POSITIONS = {
    Vector3{-1, 0, 0},   // 0
    Vector3{1, 0, 0},   // 1
    Vector3{0, 1, 0},   // 2
    Vector3{0, 0, 1},   // 3
    Vector3{0,-1, 0},   // 4
    Vector3{0, 0,-1}    // 5
};

struct TetDef {
    std::array<size_t,4> idx;
    std::array<Vector3,4> pos;
};

static const std::array<TetDef,4> TETS = {{
    {{0,1,2,3}, {RING_POSITIONS[0], RING_POSITIONS[1], RING_POSITIONS[2], RING_POSITIONS[3]}},
    {{0,1,3,4}, {RING_POSITIONS[0], RING_POSITIONS[1], RING_POSITIONS[3], RING_POSITIONS[4]}},
    {{0,1,4,5}, {RING_POSITIONS[0], RING_POSITIONS[1], RING_POSITIONS[4], RING_POSITIONS[5]}},
    {{0,1,5,2}, {RING_POSITIONS[0], RING_POSITIONS[1], RING_POSITIONS[5], RING_POSITIONS[2]}}
}};

// ── sliders ──────────────────────────────────────────────────────────────────
// Grouped: shared edge, hub-0, hub-1, ring
int ne_shared = 2;                           // {0,1}
int ne_hub0[4] = {1, 1, 1, 1};             // {0,2}, {0,3}, {0,4}, {0,5}
int ne_hub1[4] = {1, 1, 1, 1};             // {1,2}, {1,3}, {1,4}, {1,5}
int ne_ring[4] = {0, 0, 0, 0};             // {2,3}, {3,4}, {4,5}, {5,2}
int max_nc = 10;
int rng_seed = 42;
float scoop_bulge = 0.01f;
bool use_comb_merge = true;
bool use_random_ts = true;

using EdgeKey = std::pair<int,int>;

static EdgeKey make_edge(int a, int b) { return {std::min(a,b), std::max(a,b)}; }

// Build global edge counts from slider state
static std::map<EdgeKey, int> build_global_counts() {
    std::map<EdgeKey, int> counts;
    counts[{0,1}] = ne_shared;
    int ring_verts[] = {2, 3, 4, 5};
    for (int i = 0; i < 4; ++i) {
        counts[make_edge(0, ring_verts[i])] = ne_hub0[i];
        counts[make_edge(1, ring_verts[i])] = ne_hub1[i];
    }
    counts[make_edge(2,3)] = ne_ring[0];
    counts[make_edge(3,4)] = ne_ring[1];
    counts[make_edge(4,5)] = ne_ring[2];
    counts[make_edge(5,2)] = ne_ring[3];
    return counts;
}

static std::map<EdgeKey, std::vector<double>>
build_global_ts(const std::map<EdgeKey, int>& counts, int seed, bool random) {
    std::map<EdgeKey, std::vector<double>> ts;
    if (random) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> dist(0.01, 0.99);
        for (auto& [e, n] : counts) {
            std::vector<double> vals(n);
            for (int k = 0; k < n; ++k) vals[k] = dist(rng);
            std::sort(vals.begin(), vals.end());
            ts[e] = vals;
        }
    } else {
        for (auto& [e, n] : counts) {
            std::vector<double> vals(n);
            for (int k = 0; k < n; ++k) vals[k] = (k + 1.0) / (n + 1.0);
            ts[e] = vals;
        }
    }
    return ts;
}

void reconstruct_ring() {
    auto global_counts = build_global_counts();
    auto global_ts = build_global_ts(global_counts, rng_seed, use_random_ts);

    TriangleSoup::COMB_MERGE = use_comb_merge;
    TriangleSoup combined;

    for (int ti = 0; ti < 4; ++ti) {
        const auto& tet = TETS[ti];

        // Map global edges to local 6-edge frame
        std::array<int,6> local_counts{};
        std::array<std::vector<double>,6> local_ts;
        for (int ei = 0; ei < 6; ++ei) {
            int li = ALL_TET_PAIRS[ei].first;
            int lj = ALL_TET_PAIRS[ei].second;
            EdgeKey ge = make_edge(tet.idx[li], tet.idx[lj]);
            local_counts[ei] = global_counts[ge];
            local_ts[ei] = global_ts[ge];
        }

        int total = 0;
        for (int c : local_counts) total += c;
        if (total == 0) continue;

        EdgeOccupations occ;
        auto [open_c, scoop_c, normal_c] = boundary_comb_curves(
            tet.idx, local_counts, occ, false
        );
        TriangleSoup local_soup = subgrid_surface(
            tet.pos, tet.idx, local_ts,
            open_c, scoop_c, occ,
            true, (double)scoop_bulge, false
        );
        if (!local_soup.faces.empty())
            combined.add_local_soup(local_soup);
    }

    // Visualize
    polyscope::removeAllStructures();

    // Register tet wireframes
    for (int ti = 0; ti < 4; ++ti) {
        auto& tet = TETS[ti];
        std::vector<Vector3> verts(tet.pos.begin(), tet.pos.end());
        std::vector<std::vector<size_t>> faces = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
        std::string name = "Tet " + std::to_string(ti);
        polyscope::registerSurfaceMesh(name, verts, faces)
            ->setSurfaceColor({0.1f, 0.4f, 0.8f})
            ->setTransparency(0.3f)
            ->setEnabled(true);
    }

    if (!combined.faces.empty()) {
        auto* ps_mesh = polyscope::registerSurfaceMesh(
            "Ring Output", combined.vertices, combined.faces);
        ps_mesh->setSurfaceColor({0.8f, 0.2f, 0.2f});
        ps_mesh->setEnabled(true);

        // Check manifoldness
        auto [mesh, geo] = makeSurfaceMeshAndGeometry(combined.faces, combined.vertices);
        bool manifold = check_edge_manifoldness(*mesh);
        mesh->greedilyOrientFaces();
        bool oriented = mesh->isOriented();
        std::cout << "  verts: " << combined.vertices.size()
                  << "  faces: " << combined.faces.size()
                  << "  manifold: " << manifold
                  << "  oriented: " << oriented << "\n";
    } else {
        std::cout << "  (empty output)\n";
    }
}

void myCallback() {
    bool changed = false;
    changed |= ImGui::SliderInt("shared {0,1}", &ne_shared, 0, max_nc);
    changed |= ImGui::SliderInt4("hub-0 {0,2/3/4/5}", ne_hub0, 0, max_nc);
    changed |= ImGui::SliderInt4("hub-1 {1,2/3/4/5}", ne_hub1, 0, max_nc);
    changed |= ImGui::SliderInt4("ring {23/34/45/52}", ne_ring, 0, max_nc);
    changed |= ImGui::SliderInt("rng seed", &rng_seed, 0, 1000);
    changed |= ImGui::SliderFloat("scoop bulge", &scoop_bulge, 0.0f, 0.2f);
    changed |= ImGui::Checkbox("random t's", &use_random_ts);
    changed |= ImGui::Checkbox("comb merge", &use_comb_merge);
    if (changed) reconstruct_ring();
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("Interactive 4-tet ring for comb-merge debugging.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

    try { parser.ParseCLI(argc, argv); }
    catch (args::Help) { std::cout << parser; return 0; }
    catch (args::ParseError e) { std::cerr << e.what() << "\n" << parser; return 1; }

    polyscope::state::userCallback = myCallback;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::init();

    reconstruct_ring();
    polyscope::show();
    return 0;
}
