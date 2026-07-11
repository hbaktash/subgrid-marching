#include <catch2/catch_all.hpp>

#include "tet_generators.h"
#include "subgrid_MT/boundary_curve.h"

#include <set>
#include <tuple>

// ── helpers ───────────────────────────────────────────────────────────────────

// Canonical key for a CombinatorialVertex: (edge_lo, edge_hi, canonical_order).
// canonical_order normalises for the reversed-edge aliasing (matches cv::operator==).
static std::tuple<int,int,int> cv_key(const CombinatorialVertex& cv) {
    int lo = std::min(cv.i, cv.j);
    int hi = std::max(cv.i, cv.j);
    int ord = (cv.i < cv.j) ? cv.order : (cv.edge_int - cv.order + 1);
    return {lo, hi, ord};
}

static bool same_edge(const CombinatorialVertex& a, const CombinatorialVertex& b) {
    return std::min(a.i,a.j) == std::min(b.i,b.j)
        && std::max(a.i,a.j) == std::max(b.i,b.j);
}

// ── shared inputs ─────────────────────────────────────────────────────────────

static std::vector<EdgeInts> all_inputs() {
    auto v = enumerate_up_to_symmetry(10);
    auto r = random_edge_coords(30, 10000, 10);
    v.insert(v.end(), r.begin(), r.end());
    return v;
}

// ── tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("Boundary curves: every CV appears exactly once",
          "[boundary_curve]")
{
    static const std::array<size_t,4> dummy_idx = {0, 1, 2, 3};
    int tested = 0;

    for (const auto& e : all_inputs()) {
        auto counts = e.to_array();
        int total = 0;
        for (int c : counts) total += c;
        if (total == 0) continue;
        ++tested;

        INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
             <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3));

        EdgeOccupations occ;
        auto [open, scoop, normal] = boundary_comb_curves(
            dummy_idx, counts, occ, /*populate_normal=*/true
        );

        // Collect ALL CVs (open + scoop + normal) — uniqueness holds regardless
        // of even-sum: each intersection point belongs to exactly one curve.
        std::vector<std::tuple<int,int,int>> all_keys;
        for (const auto* curves : {&open, &scoop, &normal})
            for (const auto& face : *curves)
                for (const auto& cv : face.vertices)
                    all_keys.push_back(cv_key(cv));

        std::set<std::tuple<int,int,int>> key_set(all_keys.begin(), all_keys.end());

        CHECK((int)all_keys.size() == total);           // right count
        CHECK(all_keys.size() == key_set.size());       // no duplicates
    }
    REQUIRE(tested > 0);
}

TEST_CASE("Boundary curves: normal curves have no same-edge consecutive CVs",
          "[boundary_curve]")
{
    static const std::array<size_t,4> dummy_idx = {0, 1, 2, 3};
    int tested = 0;

    for (const auto& e : all_inputs()) {
        auto counts = e.to_array();
        if (std::accumulate(counts.begin(), counts.end(), 0) == 0) continue;

        EdgeOccupations occ;
        auto [open, scoop, normal] = boundary_comb_curves(
            dummy_idx, counts, occ, /*populate_normal=*/true
        );
        if (normal.empty()) continue;
        ++tested;

        INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
             <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3));

        for (const auto& curve : normal) {
            int n = (int)curve.vertices.size();
            CHECK(n >= 3);
            for (int k = 0; k < n; ++k) {
                const auto& va = curve.vertices[k];
                const auto& vb = curve.vertices[(k+1) % n];
                INFO("curve vertex " << k << ": (" << va.i << "," << va.j
                     << " ord=" << va.order << ")  next: (" << vb.i << ","
                     << vb.j << " ord=" << vb.order << ")");
                CHECK(!same_edge(va, vb));
            }
        }
    }
    REQUIRE(tested > 0);
}

TEST_CASE("Boundary curves: scoop curves have at least one same-edge consecutive pair",
          "[boundary_curve]")
{
    static const std::array<size_t,4> dummy_idx = {0, 1, 2, 3};
    int tested = 0;

    for (const auto& e : all_inputs()) {
        auto counts = e.to_array();
        if (std::accumulate(counts.begin(), counts.end(), 0) == 0) continue;

        EdgeOccupations occ;
        auto [open, scoop, normal] = boundary_comb_curves(
            dummy_idx, counts, occ, /*populate_normal=*/true
        );
        if (scoop.empty()) continue;
        ++tested;

        INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
             <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3));

        for (const auto& curve : scoop) {
            int n = (int)curve.vertices.size();
            CHECK(n >= 2);
            bool found = false;
            for (int k = 0; k < n; ++k)
                if (same_edge(curve.vertices[k], curve.vertices[(k+1) % n]))
                    { found = true; break; }
            CHECK(found);
        }
    }
    REQUIRE(tested > 0);
}
