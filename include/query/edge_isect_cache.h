#pragma once

#include "grid/grid_iterators.h"
#include "grid/explicit_tet_range.h"
#include "query/input_query_handler.h"

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

// ============================================================================
// Edge intersection reuse
// ============================================================================
//
// Each grid edge is shared by 4 (axis) or 6 (face diagonal) tets, but the
// pipeline's per-tet loop queries all six of a tet's edges from scratch, so
// every edge is queried ~5 times over the grid: 30n^3 queries against 6n^3
// unique edges. The cache serves an edge's second and later visits from memory
// instead.
//
// It is optional; with no cache the pipeline queries every tet edge directly,
// exactly as it always has.

// One edge's intersections, stored in the canonical i < j direction: t-values
// ascending in [0, 1] measured from the lower-indexed endpoint, normals
// parallel to them (empty when normals were not recorded).
struct EdgeIsect {
    std::vector<double>  ts;
    std::vector<Vector3> normals;
};

// Write a canonically-stored entry into one tet-edge slot. When the tet
// traverses the edge from its higher-indexed endpoint to its lower one, the
// t-values flip (t -> 1-t) and both arrays reverse, so the slot comes out
// ascending in the tet's own direction. See docs/construction_policy.md on edge
// parameterization. Outputs are cleared first.
void emit_edge_isect(
    const EdgeIsect& data,
    bool reversed,
    std::vector<double>& out_ts,
    std::vector<Vector3>& out_normals,
    bool record_normals
);


// ---- interface ----

class EdgeIsectCache {
public:
    virtual ~EdgeIsectCache() = default;

    // Prepare for the tets of chunk `chunk`, which for the grid is cube layer
    // iz. Chunks must be visited in increasing order; SlabEdgeCache slides its
    // window here, which is the phase that runs across threads.
    virtual void begin_chunk(
        InputQueryHandler& handler, size_t chunk,
        bool use_robust, bool record_normals
    ) = 0;

    // Fill a tet's six edge slots -- same contract as
    // InputQueryHandler::query_intersections, except that edges the cache
    // already holds are not re-queried.
    virtual void query_tet(
        InputQueryHandler& handler,
        const std::array<size_t,4>& tet_indices,
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool use_robust, bool record_normals
    ) = 0;

    // Single-edge queries actually issued to the handler. Against 30n^3 without
    // a cache, this lands at 6n^3 -- one per unique edge.
    size_t edge_queries() const { return num_edge_queries; }

protected:
    size_t num_edge_queries = 0;

    // Canonical key for an undirected edge. Node indices are (n+1)^3 < 2^32 for
    // any resolution this pipeline handles.
    static uint64_t edge_key(size_t i, size_t j) {
        const uint64_t lo = i < j ? i : j;
        const uint64_t hi = i < j ? j : i;
        return (lo << 32) | hi;
    }
};


// ---- sliding slab window (grid only, lock-free) ----
//
// Holds the edges of node planes [z, z+1] -- everything cube layer z can touch,
// and nothing else. The window's contents follow from index arithmetic alone,
// with no tet iteration and no queries, so the fill phase can pre-size its
// output and have threads write disjoint slots with no locking.
//
// Only edges that actually have intersections are stored: inside a loaded plane
// a lookup miss unambiguously means "queried, no crossings", because the plane
// was enumerated exhaustively. Memory is O(n^2), not O(n^3).
class SlabEdgeCache : public EdgeIsectCache {
public:
    SlabEdgeCache(const TetGridRange& grid_in, int num_threads_in = 1)
        : grid(grid_in), num_threads(num_threads_in) {}

    void begin_chunk(
        InputQueryHandler& handler, size_t chunk,
        bool use_robust, bool record_normals
    ) override;

    void query_tet(
        InputQueryHandler& handler,
        const std::array<size_t,4>& tet_indices,
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool use_robust, bool record_normals
    ) override;

private:
    using Map = std::unordered_map<uint64_t, EdgeIsect>;

    TetGridRange grid;
    int num_threads = 1;

    size_t cur_z = 0;      // lower node plane of the loaded window
    bool   loaded = false;

    Map lo_plane;   // both endpoints in plane cur_z
    Map cross;      // one endpoint in each plane
    Map hi_plane;   // both endpoints in plane cur_z + 1

    // Query every edge in `edges` and file the ones with hits into `into`.
    void fill(InputQueryHandler& handler,
              const std::vector<std::pair<size_t,size_t>>& edges,
              Map& into, bool use_robust, bool record_normals);

    std::vector<std::pair<size_t, size_t>> scratch;
};


// ---- building one for a tet range ----
//
// Overloaded rather than templated so each tet range says for itself whether it
// can support a cache. Returns null when `enabled` is false, and always null for
// an explicit tet mesh: the sliding window needs the grid's node planes, and
// there is nothing to gain there anyway -- those inputs arrive with their
// intersections already precomputed per edge, so PrecomputedQueryHandler is
// itself a hash lookup on the canonical edge and wrapping it in another one only
// costs time. That fallback is silent, not a warning: it is the normal outcome
// for npz input, not a misconfiguration.

std::unique_ptr<EdgeIsectCache> make_edge_cache(
    bool enabled, const TetGridRange& range, int num_threads = 1);
std::unique_ptr<EdgeIsectCache> make_edge_cache(
    bool enabled, const ExplicitTetRange& range, int num_threads = 1);

// Tets per cache chunk, or 0 when the range has no chunk structure (the cache's
// begin_chunk is then called once, for chunk 0). One cube layer for the grid:
// an edge's incident tets never span more than two consecutive layers.
inline size_t cache_chunk_size(const TetGridRange& range) { return range.tets_per_slab(); }
inline size_t cache_chunk_size(const ExplicitTetRange&)   { return 0; }
