#include <catch2/catch_all.hpp>

#include "grid/grid_iterators.h"

#include <map>
#include <set>
#include <vector>

// Structural facts about the 5-tet grid that the query cache and the parallel
// chunk driver rely on. Everything here is checked against brute-force
// iteration of the range itself, so the analytic shortcuts can never silently
// drift from what the pipeline actually visits.

namespace {

using EdgeKey = std::pair<size_t, size_t>;

EdgeKey canon(size_t i, size_t j) { return {std::min(i, j), std::max(i, j)}; }

// Every (canonical edge -> number of incident tets) pair, by walking the range.
std::map<EdgeKey, size_t> brute_force_edge_degrees(const TetGridRange& range) {
    std::map<EdgeKey, size_t> degrees;
    for (const auto& tet : range)
        for (const auto& e : ALL_TET_PAIRS)
            degrees[canon(tet.indices[e.first], tet.indices[e.second])]++;
    return degrees;
}

}  // namespace


TEST_CASE("Grid: tet_at matches sequential iteration", "[grid][random_access]") {
    // The parallel driver addresses tets by flat index instead of walking the
    // iterator; the two must agree element for element, including ordering.
    for (size_t n : {1u, 2u, 3u, 5u}) {
        TetGridRange range(n, true);

        size_t flat = 0;
        for (const auto& tet : range) {
            auto direct = range.tet_at(flat);
            INFO("n=" << n << " flat=" << flat);
            REQUIRE(direct.indices == tet.indices);
            for (int k = 0; k < 4; ++k) {
                REQUIRE(direct.positions[k].x == Catch::Approx(tet.positions[k].x));
                REQUIRE(direct.positions[k].y == Catch::Approx(tet.positions[k].y));
                REQUIRE(direct.positions[k].z == Catch::Approx(tet.positions[k].z));
            }
            ++flat;
        }
        REQUIRE(flat == range.total_tet_count());
    }
}


TEST_CASE("Grid: interior tet sits on the even-parity corners", "[grid][edges]") {
    // This is what makes the decomposition conforming: the choice of diagonal on
    // a shared cube face is a function of the lattice parity of its corners, not
    // of which cube you view it from.
    for (size_t n : {2u, 3u, 4u}) {
        TetGridRange range(n, true);
        for (size_t iz = 0; iz < n; ++iz)
        for (size_t iy = 0; iy < n; ++iy)
        for (size_t ix = 0; ix < n; ++ix) {
            auto corners = grid_cube_corner_indices(n, ix, iy, iz);
            auto tets = grid_cube_tets(ix, iy, iz, corners);
            // tets[0] is the interior (regular) tet by construction.
            for (size_t v : tets[0]) {
                auto c = grid_node_index_to_coords(n, v);
                INFO("n=" << n << " cube=(" << ix << "," << iy << "," << iz << ")");
                REQUIRE((c[0] + c[1] + c[2]) % 2 == 0);
            }
        }
    }
}


TEST_CASE("Grid: classify_edge matches brute-force degrees", "[grid][edges]") {
    for (size_t n : {1u, 2u, 3u, 4u}) {
        TetGridRange range(n, true);
        auto observed = brute_force_edge_degrees(range);

        // Every edge the grid actually uses must be classified, with the right degree.
        for (const auto& [edge, degree] : observed) {
            auto info = range.classify_edge(edge.first, edge.second);
            INFO("n=" << n << " edge=(" << edge.first << "," << edge.second << ")");
            REQUIRE(info.valid);
            REQUIRE(info.degree == degree);
        }

        // ...and nothing else may classify as valid: sweep every node pair.
        size_t valid_count = 0;
        for (size_t i = 0; i < range.node_count(); ++i)
        for (size_t j = i + 1; j < range.node_count(); ++j) {
            if (!range.classify_edge(i, j).valid) continue;
            valid_count++;
            INFO("n=" << n << " spurious edge=(" << i << "," << j << ")");
            REQUIRE(observed.count({i, j}) == 1);
        }
        REQUIRE(valid_count == observed.size());
    }
}


TEST_CASE("Grid: interior degrees are 4 (axis) and 6 (face diagonal)", "[grid][edges]") {
    // The numbers the cache's memory story is built on. n=4 is the smallest
    // resolution with a node strictly interior in every direction.
    const size_t n = 4;
    TetGridRange range(n, true);

    // A node well away from the boundary, and its 9 canonical neighbours.
    const size_t i = grid_node_coords_to_index(n, 2, 2, 2);
    size_t axis_seen = 0, diag_seen = 0;

    for (int dz = -1; dz <= 1; ++dz)
    for (int dy = -1; dy <= 1; ++dy)
    for (int dx = -1; dx <= 1; ++dx) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        const size_t j = grid_node_coords_to_index(n, 2 + dx, 2 + dy, 2 + dz);
        auto info = range.classify_edge(i, j);
        if (!info.valid) continue;
        INFO("delta=(" << dx << "," << dy << "," << dz << ")");
        if (info.isAxis) { REQUIRE(info.degree == 4); axis_seen++; }
        else             { REQUIRE(info.degree == 6); diag_seen++; }
    }

    REQUIRE(axis_seen == 6);   // +-x, +-y, +-z
    REQUIRE(diag_seen == 12);  // the even-parity diagonal of each of 12 incident faces
}


TEST_CASE("Grid: slab window is exactly the edge set of its cube layer", "[grid][edges][slab]") {
    // The cache's core invariant: the window enumerated from index arithmetic
    // covers every edge the layer's tets will look up (no misses), and adds
    // none they will not (no wasted queries).
    for (size_t n : {1u, 2u, 3u, 4u}) {
        TetGridRange range(n, true);

        for (size_t iz = 0; iz < n; ++iz) {
            std::set<EdgeKey> used;
            const size_t lo = iz * range.tets_per_slab();
            const size_t hi = lo + range.tets_per_slab();
            for (size_t flat = lo; flat < hi; ++flat) {
                auto tet = range.tet_at(flat);
                for (const auto& e : ALL_TET_PAIRS)
                    used.insert(canon(tet.indices[e.first], tet.indices[e.second]));
            }

            std::vector<EdgeKey> window;
            range.slab_window_edges(iz, window);

            INFO("n=" << n << " iz=" << iz);
            // emitted exactly once each
            std::set<EdgeKey> window_set(window.begin(), window.end());
            REQUIRE(window_set.size() == window.size());
            REQUIRE(window_set == used);
        }
    }
}


TEST_CASE("Grid: slab windows cover every edge, planes shared between neighbours",
          "[grid][edges][slab]") {
    // Together the windows must cover the whole grid (nothing is never queried),
    // and consecutive windows must overlap only in the shared node plane -- that
    // overlap is what lets the sliding window query each edge exactly once.
    const size_t n = 4;
    TetGridRange range(n, true);
    auto all_edges = brute_force_edge_degrees(range);

    std::set<EdgeKey> covered;
    std::vector<EdgeKey> prev, cur;
    for (size_t iz = 0; iz < n; ++iz) {
        range.slab_window_edges(iz, cur);
        covered.insert(cur.begin(), cur.end());

        if (iz > 0) {
            std::set<EdgeKey> prev_set(prev.begin(), prev.end()), cur_set(cur.begin(), cur.end());
            std::vector<EdgeKey> shared;
            std::set_intersection(prev_set.begin(), prev_set.end(),
                                  cur_set.begin(), cur_set.end(), std::back_inserter(shared));
            // The overlap is exactly the in-plane edges of node plane iz.
            REQUIRE(!shared.empty());
            for (const auto& e : shared) {
                auto a = grid_node_index_to_coords(n, e.first);
                auto b = grid_node_index_to_coords(n, e.second);
                INFO("iz=" << iz << " shared edge crosses planes");
                REQUIRE(a[2] == (int)iz);
                REQUIRE(b[2] == (int)iz);
            }
        }
        prev = cur;
    }

    REQUIRE(covered.size() == all_edges.size());
}
