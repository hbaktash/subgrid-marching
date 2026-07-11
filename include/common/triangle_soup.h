#pragma once

#include "nc/normalCoord_containers.h"
#include <unordered_map>

using std::vector;
using std::array;
using std::pair;


// ---- Combinatorial vertex signature ----

enum class CombVertexSigType { NORMAL, SCOOP_FACE_STEINER, SCOOP_INTERIOR_STEINER, STANDALONE };

struct CombVertexKey {
    int i = -1, j = -1;     // tet-local or global node ids, normalized i < j
    int order = -1;          // 1-based; for Steiner: min_order of the segment
    int edge_int = 0;        // total intersections on edge {i,j}; NOT part of identity, used for order flip during remap
    int k = -1;              // face vertex for SCOOP_FACE_STEINER; -1 otherwise
    CombVertexSigType type = CombVertexSigType::STANDALONE;

    // edge_int deliberately excluded from equality and hash
    bool operator==(const CombVertexKey& o) const {
        return i == o.i && j == o.j && order == o.order && k == o.k && type == o.type;
    }
};

struct CombVertexKeyHash {
    size_t operator()(const CombVertexKey& k) const {
        size_t h = std::hash<int>{}(k.i);
        h ^= std::hash<int>{}(k.j)                      + 0x9e3779b9 + (h<<6) + (h>>2);
        h ^= std::hash<int>{}(k.order)                  + 0x9e3779b9 + (h<<6) + (h>>2);
        h ^= std::hash<int>{}(k.k)                      + 0x9e3779b9 + (h<<6) + (h>>2);
        h ^= std::hash<int>{}(static_cast<int>(k.type)) + 0x9e3779b9 + (h<<6) + (h>>2);
        return h;
    }
};


// Builds a normalized (i < j) CombVertexKey from a CombinatorialVertex.
CombVertexKey key_from_cv(const CombinatorialVertex &cv);


// ---- Triangle soup ----

struct TriangleSoup {
    // Set once in main_subgrid before the pipeline loop; read-only during construction.
    static bool COMB_MERGE;

    vector<Vector3> vertices;
    vector<vector<size_t>> faces;
    vector<vector<array<size_t, 3>>> faces_per_edge;  // parallel to faces; tet face triplet per polygon edge (dual pipeline)
    vector<Vector3> dual_positions;                 // parallel to faces; one dual point per polygon (dual pipeline)
    std::unordered_map<CombVertexKey, size_t, CombVertexKeyHash> sig_to_idx;
    vector<CombVertexKey> signatures;  // parallel to vertices; populated when COMB_MERGE is true

    void add_local_soup(const vector<vector<size_t>> &local_faces, const vector<Vector3> &local_vertices);
    // Dispatches to add_local_soup_comb when COMB_MERGE is true and other has signatures; else plain concat.
    void add_local_soup(const TriangleSoup &other);
    void add_local_soup_comb(const TriangleSoup &other);

    // Remap tet-local vertex indices in all signatures: new_i = idx_map[old_i].
    // Re-normalizes i < j and flips order accordingly using edge_int.
    void remap_vertex_indices(const array<int, 4> &idx_map);

    // Remap local filtered orders to global-within-tet orders.
    // count_occupied=false: local order counts unoccupied positions (normal surface).
    // count_occupied=true:  local order counts occupied positions (evensum curve).
    void remap_orders(
        const array<vector<double>, 6> &edge_isect_ts,
        const EdgeOccupations &occ,
        bool count_occupied
    );

    // Mark all signatures whose edge touches vertex index `v` as STANDALONE.
    // Use after recursive sub-tet merge to prevent the interior centroid vertex
    // from being remapped or merged at higher levels.
    void make_vertex_standalone(int v);

    // Mark all signatures of a given type as STANDALONE.
    // Use to prevent a type from participating in further merging.
    void make_type_standalone(CombVertexSigType target_type);

    pair<vector<vector<size_t>>, vector<Vector3>>
    get_soup() const { return {faces, vertices}; }
};
