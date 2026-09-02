#pragma once

#include "common/utils.h"

#include <cstdlib>
#include <utility>

// ============================================================================
// Lattice helpers for the 5-tet grid
// ============================================================================
//
// The grid is n x n x n cubes over (n+1)^3 nodes, each cube split into 5 tets:
// one interior tet plus 4 corner tets. dx = 1/(n-1) so nodes 0..n-1 span [-1,+1]
// exactly; node n extends slightly past +1. This avoids grid-surface coincidence
// for axis-aligned SDFs.
//
// These are free functions rather than members so that the sequential iterator
// and the random-access `TetGridRange::tet_at` share one definition -- a
// divergence between the two would be a silent correctness bug.

inline double grid_dx(size_t n) {
    return n > 1 ? 1.0 / static_cast<double>(n - 1) : 1.0;
}

inline Vector3 grid_node_position(size_t n, bool center_and_normalize, int x, int y, int z) {
    const double dx = grid_dx(n);
    Vector3 p{x * dx, y * dx, z * dx};
    if (center_and_normalize) {
        Vector3 center = Vector3::constant(0.5);
        p = (p - center) * 2.0;  // fit in [-1,1]^3
    }
    return p;
}

// Node index <-> lattice coordinates. idx = z*(n+1)^2 + y*(n+1) + x.
inline std::array<int, 3> grid_node_index_to_coords(size_t n, size_t idx) {
    int z = idx / ((n + 1) * (n + 1));
    idx %= ((n + 1) * (n + 1));
    int y = idx / (n + 1);
    int x = idx % (n + 1);
    return {x, y, z};
}

inline size_t grid_node_coords_to_index(size_t n, int x, int y, int z) {
    return (size_t)z * (n + 1) * (n + 1) + (size_t)y * (n + 1) + (size_t)x;
}

// The 8 corner indices of cube (ix, iy, iz), with the following layout:
//              B----C
//             /|   /|
//            A----D |
//            | F--|-G
//            |/   |/
//            E----H
inline std::array<size_t, 8> grid_cube_corner_indices(size_t n, size_t ix, size_t iy, size_t iz) {
    return {{
        iz * (n+1) * (n+1) + iy * (n+1) + ix,              // A
        iz * (n+1) * (n+1) + iy * (n+1) + (ix+1),          // B
        iz * (n+1) * (n+1) + (iy+1) * (n+1) + (ix+1),      // C
        iz * (n+1) * (n+1) + (iy+1) * (n+1) + ix,          // D
        (iz+1) * (n+1) * (n+1) + iy * (n+1) + ix,          // E
        (iz+1) * (n+1) * (n+1) + iy * (n+1) + (ix+1),      // F
        (iz+1) * (n+1) * (n+1) + (iy+1) * (n+1) + ix,      // G
        (iz+1) * (n+1) * (n+1) + (iy+1) * (n+1) + (ix+1)   // H
    }};
}

// The 5 tets of cube (ix, iy, iz). Which of the two inscribed regular tets is
// the interior one flips with cube parity, so that the face diagonals agree
// across every shared cube face (see `TetGridRange::classify_edge`).
inline std::array<std::array<size_t, 4>, 5>
grid_cube_tets(size_t ix, size_t iy, size_t iz, const std::array<size_t, 8>& c) {
    const size_t A = c[0], B = c[1], C = c[2], D = c[3], E = c[4], F = c[5], G = c[6], H = c[7];

    if ((ix + iy + iz) % 2 == 0) {
        return {{
            {A, C, F, G},  // ACFG
            {D, A, G, C},  // DACG
            {B, A, C, F},  // BACF
            {E, F, G, A},  // EFGA
            {H, F, C, G}   // HFCG
        }};
    } else {
        return {{
            {B, E, H, D},  // BEHD
            {B, E, A, D},  // BEAD
            {B, C, H, D},  // BCHD
            {B, E, H, F},  // BEHF
            {G, E, H, D}   // GEHD
        }};
    }
}


// Lazy tet grid iterator: n³ cubes × 5 tets each, computed on-the-fly.
class TetGridIterator {
private:
    size_t n;  // grid resolution; creates n x n x n cubes
    bool center_and_normalize;

    // Current iteration state
    size_t ix, iy, iz;
    size_t tet_index;  // which tet within current cube (0-4)
    bool is_end;

public:
    // Data structure returned by dereference
    struct TetData {
        std::array<size_t, 4> indices;
        std::array<Vector3, 4> positions;
    };

    TetGridIterator(size_t n_in, bool center_normalize, bool end_iter = false)
        : n(n_in),
          center_and_normalize(center_normalize),
          ix(0), iy(0), iz(0), tet_index(0), is_end(end_iter)
    {}

    // Dereference operator: compute current tet data on-the-fly
    TetData operator*() const {
        auto corners = grid_cube_corner_indices(n, ix, iy, iz);
        auto tet_indices = grid_cube_tets(ix, iy, iz, corners)[tet_index];
        std::array<Vector3, 4> tet_positions;

        for (int i = 0; i < 4; ++i) {
            auto [x, y, z] = grid_node_index_to_coords(n, tet_indices[i]);
            tet_positions[i] = grid_node_position(n, center_and_normalize, x, y, z);
        }

        return {tet_indices, tet_positions};
    }

    // Prefix increment
    TetGridIterator& operator++() {
        if (is_end) return *this;

        tet_index++;
        if (tet_index >= 5) {
            tet_index = 0;
            ix++;
            if (ix >= n) {
                ix = 0;
                iy++;
                if (iy >= n) {
                    iy = 0;
                    iz++;
                    if (iz >= n) {
                        is_end = true;
                    }
                }
            }
        }
        return *this;
    }

    bool operator==(const TetGridIterator& other) const {
        return is_end == other.is_end &&
               (is_end || (ix == other.ix && iy == other.iy && iz == other.iz && tet_index == other.tet_index));
    }

    bool operator!=(const TetGridIterator& other) const {
        return !(*this == other);
    }

    size_t total_tet_count() const { return static_cast<size_t>(n) * n * n * 5; }
};


// ---- Edge classification ----

// What `TetGridRange::classify_edge` reports about a node pair.
struct GridEdgeInfo {
    bool   valid  = false;  // false when {i, j} is not an edge of this grid
    bool   isAxis = false;  // true: axis-aligned edge; false: face diagonal
    size_t degree = 0;      // number of tets of the grid containing the edge
};


// Range wrapper for range-based for loops
class TetGridRange {
    size_t n;
    bool center_and_normalize;

    // How many cube layers along one axis contain the lattice coordinate
    // `coord`: cubes with index in {coord-1, coord}, clamped to [0, n-1].
    size_t cube_span(int coord) const {
        int lo = std::max(0, coord - 1);
        int hi = std::min((int)n - 1, coord);
        return hi >= lo ? (size_t)(hi - lo + 1) : 0;
    }

public:
    using TetData = TetGridIterator::TetData;

    TetGridRange(size_t n_in, bool cnorm = true)
        : n(n_in), center_and_normalize(cnorm) {}

    TetGridIterator begin() const {
        return TetGridIterator(n, center_and_normalize, false);
    }

    TetGridIterator end() const {
        return TetGridIterator(n, center_and_normalize, true);
    }

    size_t resolution() const { return n; }
    size_t node_count() const { return (n + 1) * (n + 1) * (n + 1); }

    size_t total_tet_count() const {
        return n * n * n * 5;
    }

    // ---- random access ----
    //
    // Flat tet indexing in exactly the order `operator++` advances (tet_index
    // fastest, then ix, then iy, then iz), so a contiguous flat range is a
    // contiguous run of the sequential iteration -- and one cube layer `iz` is
    // the flat range [iz * tets_per_slab(), (iz+1) * tets_per_slab()).
    size_t tets_per_slab() const { return n * n * 5; }
    size_t slab_count() const { return n; }

    TetData tet_at(size_t flat_index) const {
        const size_t tet_index = flat_index % 5;
        size_t cube = flat_index / 5;
        const size_t ix = cube % n;
        cube /= n;
        const size_t iy = cube % n;
        const size_t iz = cube / n;

        auto corners = grid_cube_corner_indices(n, ix, iy, iz);
        auto tet_indices = grid_cube_tets(ix, iy, iz, corners)[tet_index];
        std::array<Vector3, 4> tet_positions;
        for (int i = 0; i < 4; ++i) {
            auto [x, y, z] = grid_node_index_to_coords(n, tet_indices[i]);
            tet_positions[i] = grid_node_position(n, center_and_normalize, x, y, z);
        }
        return {tet_indices, tet_positions};
    }

    // ---- edge classification ----
    //
    // Every edge of the 5-tet grid is either an axis edge or a face diagonal;
    // there are no body diagonals. The interior tet of a cube sits on the four
    // corners whose lattice parity (x+y+z) is even -- that characterization is
    // global, which is why the diagonal picked on a shared cube face is the same
    // from both sides and the mesh is conforming.
    //
    // Degrees follow from two per-cube facts: each of the 12 axis edges lies in
    // exactly 1 of the 5 tets (4 corner tets x 3 axis edges = 12), and each of
    // the 6 face diagonals lies in 3 (the interior tet plus 2 corner tets).
    // Multiplying by the number of cubes sharing the edge gives 4 for an
    // interior axis edge and 6 for an interior face diagonal.
    GridEdgeInfo classify_edge(size_t i, size_t j) const {
        GridEdgeInfo info;
        if (i == j || i >= node_count() || j >= node_count()) return info;

        auto a = grid_node_index_to_coords(n, i);
        auto b = grid_node_index_to_coords(n, j);

        int steps = 0, free_axis = -1;
        for (int k = 0; k < 3; ++k) {
            int d = std::abs(b[k] - a[k]);
            if (d == 0) free_axis = k;
            else if (d == 1) steps++;
            else return info;  // stride > 1: never an edge
        }

        if (steps == 1) {
            // Axis edge: shared by every cube spanning it in the two
            // perpendicular directions, 1 tet in each.
            info.valid = true;
            info.isAxis = true;
            info.degree = 1;
            for (int k = 0; k < 3; ++k)
                if (a[k] == b[k]) info.degree *= cube_span(a[k]);
            return info;
        }

        if (steps == 2) {
            // Face diagonal. Only the one joining the two even-parity corners of
            // the face is used; the other diagonal of that face is not an edge.
            // Both endpoints share a parity here (the delta has two 1s), so one
            // check settles it.
            if (((a[0] + a[1] + a[2]) & 1) != 0) return info;
            info.valid = true;
            info.isAxis = false;
            info.degree = 3 * cube_span(a[free_axis]);
            return info;
        }

        return info;  // steps == 3: body diagonal, not an edge of this decomposition
    }

    // Number of tets of the grid containing edge {i, j}; 0 when it is not an edge.
    size_t edge_tet_degree(size_t i, size_t j) const {
        return classify_edge(i, j).degree;
    }

    // ---- slab windows ----
    //
    // Every tet of cube layer `iz` has all four corners in node planes iz and
    // iz+1, so the edges it can query are exactly the grid edges with both
    // endpoints in those two planes -- derivable from index arithmetic alone,
    // with no tet iteration and no intersection queries. That is what lets the
    // parallel query phase pre-size its output and write disjoint slots.
    //
    // Writes canonical (min, max) node pairs; each edge appears exactly once.
    //
    // The window splits into three disjoint parts -- the in-plane edges of each
    // of the two planes, and the edges crossing between them. Consecutive
    // windows share exactly one of those parts (plane iz+1 becomes plane iz), so
    // a sliding cache carries it over instead of re-querying it, and every edge
    // of the grid is queried exactly once.
    void slab_window_edges(size_t iz, std::vector<std::pair<size_t, size_t>>& out) const {
        out.clear();
        if (iz >= n) return;
        std::vector<std::pair<size_t, size_t>> part;
        plane_edges(iz, part);      out.insert(out.end(), part.begin(), part.end());
        cross_edges(iz, part);      out.insert(out.end(), part.begin(), part.end());
        plane_edges(iz + 1, part);  out.insert(out.end(), part.begin(), part.end());
    }

    // Edges with both endpoints in node plane z.
    void plane_edges(size_t z, std::vector<std::pair<size_t, size_t>>& out) const {
        out.clear();
        if (z > n) return;
        // First nonzero component positive, so each edge is emitted once.
        static constexpr int DIRS[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
        for (int y = 0; y <= (int)n; ++y)
        for (int x = 0; x <= (int)n; ++x)
            for (const auto& d : DIRS) {
                const int bx = x + d[0], by = y + d[1];
                if (bx < 0 || bx > (int)n || by < 0 || by > (int)n) continue;
                const size_t i = grid_node_coords_to_index(n, x, y, (int)z);
                const size_t j = grid_node_coords_to_index(n, bx, by, (int)z);
                if (!classify_edge(i, j).valid) continue;
                out.emplace_back(std::min(i, j), std::max(i, j));
            }
    }

    // Edges with one endpoint in node plane z and the other in plane z+1.
    void cross_edges(size_t z, std::vector<std::pair<size_t, size_t>>& out) const {
        out.clear();
        if (z >= n + 1) return;
        // Emitted from the lower plane only, so each edge is emitted once even
        // though the in-plane deltas are not sign-canonical here.
        static constexpr int DIRS[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        for (int y = 0; y <= (int)n; ++y)
        for (int x = 0; x <= (int)n; ++x)
            for (const auto& d : DIRS) {
                const int bx = x + d[0], by = y + d[1];
                if (bx < 0 || bx > (int)n || by < 0 || by > (int)n) continue;
                const size_t i = grid_node_coords_to_index(n, x, y, (int)z);
                const size_t j = grid_node_coords_to_index(n, bx, by, (int)z + 1);
                if (!classify_edge(i, j).valid) continue;
                out.emplace_back(std::min(i, j), std::max(i, j));
            }
    }
};
