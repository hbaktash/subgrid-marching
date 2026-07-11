#include "query/input_query_handler.h"
#include "query/intersection_query.h"
#include "query/sdf_queries.h"
#include "assembly/mesh_processing.h"
#include "sdf/sdf.hpp"

#include <iostream>
#include <random>

// ============================================================================
// SDF Handler Implementation
// ============================================================================



// query count logging utilities
void InputQueryHandler::update_global_query_count_map(
    std::array<size_t, 4>& global_tet_indices,
    bool active_cell
){
    // Update global edge query counts
    for (size_t edge_idx = 0; edge_idx < 6; ++edge_idx){
        if (local_edge_query_counts[edge_idx] == 0) continue; // skip if no queries for this edge
        int i = ALL_TET_PAIRS[edge_idx].first;
        int j = ALL_TET_PAIRS[edge_idx].second;
        std::pair<int,int> edge_key = std::make_pair(
            std::min(global_tet_indices[i], global_tet_indices[j]),
            std::max(global_tet_indices[i], global_tet_indices[j])
        );
        if (active_cell) {
            active_edge_query_counts[edge_key] = local_edge_query_counts[edge_idx];
        }
        edge_query_counts[edge_key] = local_edge_query_counts[edge_idx];
    }
    // refresh local counts for next query
    local_edge_query_counts = {0, 0, 0, 0, 0, 0};
}

SDFQueryHandler::SDFQueryHandler(const std::string& name, float step_size)
    : sdf_name(name), min_step_size(step_size)
{
    // Create the SDF evaluation function
    sdf_func = [name](const Vector3& p) -> float {
        glm::vec3 glm_p(p.x, p.y, p.z);
        return sdf::evaluate(name, glm_p);
    };
    
    std::cout << "Loaded SDF from file: " << sdf_name << std::endl;
}

void SDFQueryHandler::query_intersections(
    const std::array<Vector3,4>& tet_positions,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool /*useRobust*/,
    bool /*recordNormals*/
) {
    std::array<size_t, 6> query_count_per_edge;
    tet_edge_intersections_SDF(
        tet_positions, sdf_func,
        edge_isect_ts, edge_isect_normals,
        query_count_per_edge,
        min_step_size
    );
    total_queries += query_count_per_edge[0] + query_count_per_edge[1] + query_count_per_edge[2] + query_count_per_edge[3] + query_count_per_edge[4] + query_count_per_edge[5];
    // Update edge query counts for logging
    local_edge_query_counts = query_count_per_edge;
}



void SDFQueryHandler::query_normal(
    const Vector3& q,
    Vector3& normal,
    bool verbose
) {
    normal = query_normal_SDF(q, sdf_func);
}

// ============================================================================
// Mesh Handler Implementation
// ============================================================================

MeshQueryHandler::MeshQueryHandler(
    const std::vector<Vector3>& pos,
    const std::vector<std::vector<size_t>>& polys
) : positions(pos), polygons(polys)
{
    // Build accelerator structure once
    accel = get_fcpw_accel(polygons, positions);
    std::cout << "Built FCPW accelerator for mesh with " << positions.size() 
              << " vertices and " << polygons.size() << " faces." << std::endl;
}

void MeshQueryHandler::query_intersections(
    const std::array<Vector3,4>& tet_positions,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool useRobust,
    bool recordNormals
) {
    find_single_tet_edge_intersections_fcpw(
        tet_positions, accel,
        edge_isect_ts, edge_isect_normals,
        useRobust, recordNormals
    );
}

void MeshQueryHandler::query_normal(
    const Vector3& q,
    Vector3& normal,
    bool verbose
) {
    Vector3 closest_point;
    closest_point_fcpw(q, accel, closest_point, normal);
}

// ============================================================================
// Preprocessing Utilities
// ============================================================================

PreprocessedMeshData preprocess_input_mesh(
    SimplePolygonMesh& mesh,
    double normalize_scale
) {
    // Weld coincident vertices first so that soup representations of watertight
    // meshes stay watertight (duplicated corners collapse to one vertex).
    mesh.mergeIdenticalVertices();
    mesh.triangulate();

    PreprocessedMeshData result;
    result.positions = mesh.vertexCoordinates;
    result.polygons = mesh.polygons;

    // Center and normalize the mesh to fit inside the grid box.
    center_and_normalize(result.positions, normalize_scale);

    // Decorrelate the mesh from the axis-aligned grid with a single fixed-seed
    // rigid translation. This replaces per-vertex jitter: the mesh stays
    // geometrically exact (so welding/watertightness is preserved) while grid
    // planes no longer coincide with axis-aligned mesh faces. Applied AFTER
    // normalization, otherwise the re-centering would cancel it. The magnitude
    // is half the grid-boundary clearance, so the mesh (max |coord| <=
    // normalize_scale) is guaranteed to stay strictly inside the grid.
    {
        const double t_mag = (1.0 - normalize_scale) * 0.5; // 0.01 at scale 0.98
        std::mt19937 rng(1u);                                // fixed seed for reproducibility
        std::normal_distribution<double> gaussian(0.0, 1.0);
        Vector3 dir{gaussian(rng), gaussian(rng), gaussian(rng)};
        double nrm = dir.norm();
        Vector3 t = (nrm > 1e-12 ? dir / nrm : Vector3{1.0, 0.0, 0.0}) * t_mag;
        for (Vector3& p : result.positions) p += t;
    }

    return result;
}
