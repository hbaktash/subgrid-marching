#pragma once

#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "fcpw/fcpw.h"

#include <vector>
#include <array>
#include <functional>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ============================================================================
// Base Interface
// ============================================================================

class InputQueryHandler {
public:
    virtual ~InputQueryHandler() = default;

    // Diagnostic counters. These are plain (non-atomic) shared state, so a
    // handler driven from several threads must have collect_query_stats off;
    // the parallel pipeline clears it and skips update_global_query_count_map.
    bool collect_query_stats = true;
    size_t total_queries = 0; // for logging and debugging
    // map from edges to query counts for detailed logging
    std::array<size_t, 6> local_edge_query_counts = {0, 0, 0, 0, 0, 0}; // for current query
    std::map<std::pair<int,int>, size_t> edge_query_counts;
    std::map<std::pair<int,int>, size_t> active_edge_query_counts;
    
    // subgrid-style intersection queries: returns edge intersection t-values and normals
    virtual void query_intersections(
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) = 0;

    // Index-aware overload: default delegates to position-only version.
    // PrecomputedQueryHandler overrides this to use indices for lookup.
    virtual void query_intersections(
        const std::array<size_t,4>& tet_indices,
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    );

    // Ask for every edge in the canonical min(i,j) -> max(i,j) direction and flip
    // the result for tets that traverse it the other way, instead of querying in
    // each tet's own local direction.
    //
    // It matters because no handler is exactly direction-symmetric: FCPW works in
    // float internally, and the SDF march samples at different points depending
    // on which end it starts from. Without this, two tets sharing a grid edge can
    // disagree about where -- and for SDF input, how many times -- the surface
    // crosses it. With it they always agree, and a cached run matches an uncached
    // one bit for bit.
    //
    // Costs no extra queries, only a reversal on the edges that need flipping.
    // The pipelines set this from their `canonical_queries` option, which is on
    // by default; it is off here so a handler used directly keeps the raw
    // per-call behaviour.
    bool canonical_edge_queries = false;

    // Single-edge intersection query along pi -> pj: t-values ascending in
    // [0, 1] measured from pi, normals parallel to them. Outputs are cleared
    // first. `query_intersections` is six of these; splitting them out is what
    // lets the edge cache (query/edge_isect_cache.h) query each grid edge once
    // instead of once per incident tet. Handlers pass whichever of the index /
    // position pair they actually key on.
    virtual bool supports_edge_query() const { return false; }
    virtual void query_edge(
        size_t global_i,
        size_t global_j,
        const Vector3& pi,
        const Vector3& pj,
        std::vector<double>& out_ts,
        std::vector<Vector3>& out_normals,
        bool useRobust = false,
        bool recordNormals = true
    );

    // Single-point normal query
    virtual void query_normal(
        const Vector3& q, 
        Vector3& normal, 
        bool verbose = false
    ) = 0;

    // Metadata
    virtual bool is_sdf() const = 0;
    virtual bool is_mesh() const = 0;
    
    // Visualization support - for registering with polyscope
    virtual bool has_mesh_data() const = 0;
    virtual const std::vector<Vector3>& get_mesh_positions() const = 0;
    virtual const std::vector<std::vector<size_t>>& get_mesh_polygons() const = 0;

    // query count logging utilities
    void update_global_query_count_map(
        std::array<size_t, 4>& global_tet_indices, 
        bool active_cell
    );
};

// Scoped override of a handler's diagnostic/canonicalization mode, restored on
// exit so a handler reused across pipeline runs is left as it was found.
struct QueryModeGuard {
    InputQueryHandler& h;
    bool prev_canonical, prev_stats;
    QueryModeGuard(InputQueryHandler& handler, bool canonical, bool stats)
        : h(handler), prev_canonical(handler.canonical_edge_queries),
          prev_stats(handler.collect_query_stats) {
        h.canonical_edge_queries = canonical;
        h.collect_query_stats = stats;
    }
    ~QueryModeGuard() {
        h.canonical_edge_queries = prev_canonical;
        h.collect_query_stats = prev_stats;
    }
};

// ============================================================================
// SDF Handler
// ============================================================================

class SDFQueryHandler : public InputQueryHandler {
private:
    std::string sdf_name;
    std::function<float(const Vector3&)> sdf_func;
    float min_step_size;
    
public:
    SDFQueryHandler(const std::string& name, float step_size = 1e-5f);

    // Names of all built-in SDFs (sorted), and a membership test. Backed by the
    // sdf-dataset registry, so they never go out of sync with what evaluate() accepts.
    static std::vector<std::string> available_sdf_names();
    static bool is_valid_sdf_name(const std::string& name);
    
    void query_intersections(
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;
    
    bool supports_edge_query() const override { return true; }
    void query_edge(
        size_t global_i, size_t global_j,
        const Vector3& pi, const Vector3& pj,
        std::vector<double>& out_ts,
        std::vector<Vector3>& out_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;

    void query_normal(
        const Vector3& q,
        Vector3& normal,
        bool verbose = false
    ) override;

    bool is_sdf() const override { return true; }
    bool is_mesh() const override { return false; }
    bool has_mesh_data() const override { return false; }
    
    const std::vector<Vector3>& get_mesh_positions() const override {
        throw std::runtime_error("SDFQueryHandler: no mesh data available");
    }
    const std::vector<std::vector<size_t>>& get_mesh_polygons() const override {
        throw std::runtime_error("SDFQueryHandler: no mesh data available");
    }

};

// ============================================================================
// Mesh Handler
// ============================================================================

class MeshQueryHandler : public InputQueryHandler {
private:
    std::vector<Vector3> positions;
    std::vector<std::vector<size_t>> polygons;
    fcpw::Scene<3> accel;
    
public:
    // Constructor takes preprocessed mesh data and builds accelerator
    MeshQueryHandler(
        const std::vector<Vector3>& pos,
        const std::vector<std::vector<size_t>>& polys
    );
    
    void query_intersections(
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;
    
    bool supports_edge_query() const override { return true; }
    void query_edge(
        size_t global_i, size_t global_j,
        const Vector3& pi, const Vector3& pj,
        std::vector<double>& out_ts,
        std::vector<Vector3>& out_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;

    void query_normal(
        const Vector3& q,
        Vector3& normal,
        bool verbose = false
    ) override;

    bool is_sdf() const override { return false; }
    bool is_mesh() const override { return true; }
    bool has_mesh_data() const override { return true; }

    const std::vector<Vector3>& get_mesh_positions() const override { return positions; }
    const std::vector<std::vector<size_t>>& get_mesh_polygons() const override { return polygons; }
};

// ============================================================================
// Preprocessing Utilities
// ============================================================================

struct PreprocessedMeshData {
    std::vector<Vector3> positions;
    std::vector<std::vector<size_t>> polygons;
};

// Preprocess a raw input mesh in place: weld coincident vertices, triangulate,
// normalize to fit inside the grid box, then apply a single rigid translation
// (seeded by `seed`) to decorrelate the mesh from the axis-aligned grid. No
// per-vertex jitter, so the mesh stays geometrically exact (welding/watertightness
// kept). `seed` controls the translation direction only; the default reproduces
// the historical fixed offset.
PreprocessedMeshData preprocess_input_mesh(
    SimplePolygonMesh& mesh,
    double normalize_scale = 0.98,
    unsigned int seed = 1u
);
