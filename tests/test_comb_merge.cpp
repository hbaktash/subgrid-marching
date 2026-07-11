#include <catch2/catch_all.hpp>

#include "tet_generators.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/primal_reconstruction.h"
#include "assembly/mesh_processing.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include <numeric>

// ── helpers ──────────────────────────────────────────────────────────────────

struct CombMergeGuard {
    bool prev;
    CombMergeGuard(bool val) : prev(TriangleSoup::COMB_MERGE) { TriangleSoup::COMB_MERGE = val; }
    ~CombMergeGuard() { TriangleSoup::COMB_MERGE = prev; }
};

static std::vector<EdgeInts> all_inputs() {
    auto v = enumerate_up_to_symmetry(10);
    auto r = random_edge_coords(30, 10000, 10, 42);
    v.insert(v.end(), r.begin(), r.end());
    return v;
}

// ── tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("Single-tet: local soup is manifold",
          "[comb_merge][single_tet]")
{
    static const auto tet_pos = standard_tet_positions();
    static const std::array<size_t,4> idx = {0, 1, 2, 3};
    std::mt19937 rng(123);
    int tested = 0;

    for (const auto& e : all_inputs()) {
        auto counts = e.to_array();
        int total = std::accumulate(counts.begin(), counts.end(), 0);
        if (total == 0) continue;

        for (int trial = 0; trial < 4; ++trial) {
            auto ts = random_isect_ts(e, rng);

            EdgeOccupations occ;
            auto [open, scoop, normal] = boundary_comb_curves(idx, counts, occ, false);

            TriangleSoup soup = subgrid_surface(
                tet_pos, idx, ts,
                open, scoop, occ,
                true, 0.01, false
            );
            if (soup.faces.empty()) continue;
            ++tested;

            INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
                 <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3)
                 << "  trial: " << trial);

            auto [mesh, geo] = makeSurfaceMeshAndGeometry(soup.faces, soup.vertices);
            REQUIRE(mesh != nullptr);
            CHECK(check_edge_manifoldness(*mesh));

            if (is_even_sum(e)) {
                mesh->greedilyOrientFaces();
                CHECK(mesh->isOriented());
            }
        }
    }
    REQUIRE(tested > 0);
}

TEST_CASE("Two-tet comb-merge: shared face produces manifold soup",
          "[comb_merge][two_tet]")
{
    // Tet A = {0,1,2,3}, Tet B = {0,1,2,4} — share face {0,1,2}
    // Shared edges: {0,1}, {0,2}, {1,2} (indices 0,1,3 in ALL_TET_PAIRS)
    static const std::array<Vector3,5> all_pos = {
        Vector3{0,0,0}, Vector3{1,0,0}, Vector3{0,1,0}, Vector3{0,0,1}, Vector3{0,0,-1}
    };
    static const std::array<size_t,4> idx_A = {0, 1, 2, 3};
    static const std::array<size_t,4> idx_B = {0, 1, 2, 4};
    static const std::array<Vector3,4> pos_A = {all_pos[0], all_pos[1], all_pos[2], all_pos[3]};
    static const std::array<Vector3,4> pos_B = {all_pos[0], all_pos[1], all_pos[2], all_pos[4]};

    CombMergeGuard guard(true);
    std::mt19937 rng(77);
    int tested = 0;

    for (const auto& e_A : all_inputs()) {
        auto counts_A = e_A.to_array();
        int total_A = std::accumulate(counts_A.begin(), counts_A.end(), 0);
        if (total_A == 0) continue;

        for (int trial = 0; trial < 4; ++trial) {
            // Generate t-values for tet A
            auto ts_A = random_isect_ts(e_A, rng);

            // Build edge ints for tet B: shared edges same count, non-shared random
            std::uniform_int_distribution<int> dist(0, 4);
            std::array<int,6> counts_B_arr{};
            counts_B_arr[0] = counts_A[0]; // edge {0,1} shared
            counts_B_arr[1] = counts_A[1]; // edge {0,2} shared
            counts_B_arr[2] = dist(rng);   // edge {0,4} independent
            counts_B_arr[3] = counts_A[3]; // edge {1,2} shared
            counts_B_arr[4] = dist(rng);   // edge {1,4} independent
            counts_B_arr[5] = dist(rng);   // edge {2,4} independent
            EdgeInts e_B(counts_B_arr);

            // Generate t-values for tet B: shared edges use A's values
            auto ts_B = random_isect_ts(e_B, rng);
            ts_B[0] = ts_A[0]; // edge {0,1}
            ts_B[1] = ts_A[1]; // edge {0,2}
            ts_B[3] = ts_A[3]; // edge {1,2}

            // Reconstruct tet A
            EdgeOccupations occ_A;
            auto [open_A, scoop_A, normal_A] = boundary_comb_curves(idx_A, counts_A, occ_A, false);
            TriangleSoup soup_A = subgrid_surface(pos_A, idx_A, ts_A, open_A, scoop_A, occ_A, true, 0.01, false);

            // Reconstruct tet B
            EdgeOccupations occ_B;
            auto [open_B, scoop_B, normal_B] = boundary_comb_curves(idx_B, counts_B_arr, occ_B, false);
            TriangleSoup soup_B = subgrid_surface(pos_B, idx_B, ts_B, open_B, scoop_B, occ_B, true, 0.01, false);

            if (soup_A.faces.empty() && soup_B.faces.empty()) continue;
            ++tested;

            // Comb-merge both into one soup
            TriangleSoup combined;
            combined.add_local_soup(soup_A);
            combined.add_local_soup(soup_B);

            if (combined.faces.empty()) continue;

            INFO("EdgeInts A: " << e_A.get(0,1) <<" "<< e_A.get(0,2) <<" "<< e_A.get(0,3)
                 <<" "<< e_A.get(1,2) <<" "<< e_A.get(1,3) <<" "<< e_A.get(2,3)
                 << "  trial: " << trial);
            INFO("EdgeInts B: " << e_B.get(0,1) <<" "<< e_B.get(0,2) <<" "<< e_B.get(0,3)
                 <<" "<< e_B.get(1,2) <<" "<< e_B.get(1,3) <<" "<< e_B.get(2,3));

            auto [mesh, geo] = makeSurfaceMeshAndGeometry(combined.faces, combined.vertices);
            REQUIRE(mesh != nullptr);
            CHECK(check_edge_manifoldness(*mesh));
        }
    }
    REQUIRE(tested > 0);
}

TEST_CASE("Four-tet comb-merge: edge ring produces manifold soup",
          "[comb_merge][four_tet]")
{
    // 4 tets sharing edge {0,1}, arranged in a ring:
    //   tet 0: {0,1,2,3}   tet 1: {0,1,3,4}   tet 2: {0,1,4,5}   tet 3: {0,1,5,2}
    // Vertices 2,3,4,5 form a ring around edge {0,1}.
    // Global edges (13 total):
    //   shared by all 4: {0,1}
    //   shared by 2 tets: {0,2},{1,2} (tets 0&3), {0,3},{1,3} (tets 0&1),
    //                     {0,4},{1,4} (tets 1&2), {0,5},{1,5} (tets 2&3)
    //   internal to 1 tet: {2,3} (tet 0), {3,4} (tet 1), {4,5} (tet 2), {5,2} (tet 3)

    using Edge = std::pair<int,int>;

    static const std::array<Vector3,6> all_pos = {
        Vector3{0, 0, 0},    // 0
        Vector3{1, 0, 0},    // 1
        Vector3{0, 1, 0},    // 2
        Vector3{0, 0, 1},    // 3
        Vector3{0,-1, 0},    // 4
        Vector3{0, 0,-1}     // 5
    };

    struct TetDef {
        std::array<size_t,4> idx;
        std::array<Vector3,4> pos;
    };
    static const std::array<TetDef,4> tets = {{
        {{0,1,2,3}, {all_pos[0], all_pos[1], all_pos[2], all_pos[3]}},
        {{0,1,3,4}, {all_pos[0], all_pos[1], all_pos[3], all_pos[4]}},
        {{0,1,4,5}, {all_pos[0], all_pos[1], all_pos[4], all_pos[5]}},
        {{0,1,5,2}, {all_pos[0], all_pos[1], all_pos[5], all_pos[2]}}
    }};

    CombMergeGuard guard(true);
    std::mt19937 rng(99);
    std::uniform_int_distribution<int> count_dist(0, 20);
    int tested = 0;

    for (int trial = 0; trial < 4000; ++trial) {
        // Generate one count per global edge
        std::map<Edge, int> global_counts;
        std::map<Edge, std::vector<double>> global_ts;

        auto get_edge = [](int a, int b) -> Edge {
            return {std::min(a,b), std::max(a,b)};
        };

        // Collect all global edges from the 4 tets
        std::set<Edge> all_edges;
        for (const auto& tet : tets)
            for (int ei = 0; ei < 6; ++ei) {
                int li = ALL_TET_PAIRS[ei].first;
                int lj = ALL_TET_PAIRS[ei].second;
                all_edges.insert(get_edge(tet.idx[li], tet.idx[lj]));
            }

        // Assign random counts and t-values per global edge
        std::uniform_real_distribution<double> t_dist(0.01, 0.99);
        for (auto& e : all_edges) {
            int n = count_dist(rng);
            global_counts[e] = n;
            std::vector<double> ts(n);
            for (int k = 0; k < n; ++k) ts[k] = t_dist(rng);
            std::sort(ts.begin(), ts.end());
            global_ts[e] = ts;
        }

        // Reconstruct each tet and comb-merge into one soup
        TriangleSoup combined;
        bool any_output = false;

        for (const auto& tet : tets) {
            // Build local counts and t-values from global map
            std::array<int,6> local_counts{};
            std::array<std::vector<double>,6> local_ts;
            for (int ei = 0; ei < 6; ++ei) {
                int li = ALL_TET_PAIRS[ei].first;
                int lj = ALL_TET_PAIRS[ei].second;
                Edge ge = get_edge(tet.idx[li], tet.idx[lj]);
                local_counts[ei] = global_counts[ge];
                local_ts[ei] = global_ts[ge];
            }

            int total = std::accumulate(local_counts.begin(), local_counts.end(), 0);
            if (total == 0) continue;

            EdgeOccupations occ;
            auto [open_c, scoop_c, normal_c] = boundary_comb_curves(
                tet.idx, local_counts, occ, false
            );
            TriangleSoup local_soup = subgrid_surface(
                tet.pos, tet.idx, local_ts,
                open_c, scoop_c, occ, true, 0.01, false
            );
            if (!local_soup.faces.empty()) {
                combined.add_local_soup(local_soup);
                any_output = true;
            }
        }

        if (!any_output || combined.faces.empty()) continue;
        ++tested;

        std::string edge_info;
        for (auto& [e, c] : global_counts)
            if (c > 0) edge_info += "{" + std::to_string(e.first) + "," + std::to_string(e.second) + "}=" + std::to_string(c) + " ";
        INFO("trial: " << trial << "  edges: " << edge_info);

        auto [mesh, geo] = makeSurfaceMeshAndGeometry(combined.faces, combined.vertices);
        REQUIRE(mesh != nullptr);
        CHECK(check_edge_manifoldness(*mesh));
    }
    REQUIRE(tested > 0);
}
