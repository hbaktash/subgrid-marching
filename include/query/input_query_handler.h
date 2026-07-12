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
    ) {
        query_intersections(tet_positions, edge_isect_ts, edge_isect_normals, useRobust, recordNormals);
    }
    
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
// normalize to fit inside the grid box, then apply a single fixed-seed rigid
// translation to decorrelate the mesh from the axis-aligned grid. No per-vertex
// jitter, so the mesh stays geometrically exact (welding/watertightness kept).
PreprocessedMeshData preprocess_input_mesh(
    SimplePolygonMesh& mesh,
    double normalize_scale = 0.98
);
