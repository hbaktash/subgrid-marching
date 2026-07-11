#pragma once

#include "query/input_query_handler.h"
#include "common/npz_reader.h"
#include <unordered_map>

// Query handler that looks up precomputed edge intersections from an NPZ file.
// Edges are stored canonically as (min(i,j), max(i,j)) with t-values parameterized
// from the smaller to the larger index. When a tet's local edge runs in the reverse
// direction, t-values are flipped (1-t) and reversed.

class PrecomputedQueryHandler : public InputQueryHandler {
public:
    struct EdgeData {
        std::vector<double> ts;
        std::vector<Vector3> normals;  // empty if normals not provided
    };

private:
    struct PairHash {
        size_t operator()(const std::pair<size_t,size_t>& p) const {
            size_t h = std::hash<size_t>{}(p.first);
            h ^= std::hash<size_t>{}(p.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    std::unordered_map<std::pair<size_t,size_t>, EdgeData, PairHash> edge_map;
    bool has_normals = false;

public:
    // Build from NPZ arrays:
    //   edges (E,2) int32 — must be stored with smaller index first
    //   isect_offsets (E+1,) int32 — CSR offsets
    //   isect_ts (N,) float64 — packed t-values
    //   isect_normals (N,3) float64 — optional
    PrecomputedQueryHandler(
        const NpyArray& edges_arr,
        const NpyArray& offsets_arr,
        const NpyArray& ts_arr,
        const NpyArray* normals_arr = nullptr
    );

    void query_intersections(
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;

    void query_intersections(
        const std::array<size_t,4>& tet_indices,
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;

    void query_normal(const Vector3& q, Vector3& normal, bool verbose = false) override;

    bool is_sdf() const override { return false; }
    bool is_mesh() const override { return false; }
    bool has_mesh_data() const override { return false; }
    const std::vector<Vector3>& get_mesh_positions() const override {
        throw std::runtime_error("PrecomputedQueryHandler: no mesh data");
    }
    const std::vector<std::vector<size_t>>& get_mesh_polygons() const override {
        throw std::runtime_error("PrecomputedQueryHandler: no mesh data");
    }
};
