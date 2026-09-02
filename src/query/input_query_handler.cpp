#include "query/input_query_handler.h"
#include "query/edge_isect_cache.h"
#include "query/intersection_query.h"
#include "query/sdf_queries.h"
#include "common/utils.h"
#include "sdf/sdf.hpp"

#include <iostream>
#include <random>
#include <algorithm>

// ============================================================================
// SDF Handler Implementation
// ============================================================================



void InputQueryHandler::query_intersections(
    const std::array<size_t,4>& tet_indices,
    const std::array<Vector3,4>& tet_positions,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool useRobust,
    bool recordNormals
){
    if (!canonical_edge_queries || !supports_edge_query()) {
        query_intersections(tet_positions, edge_isect_ts, edge_isect_normals, useRobust, recordNormals);
        return;
    }

    EdgeIsect scratch;
    for (int e = 0; e < 6; ++e) {
        const int li = ALL_TET_PAIRS[e].first, lj = ALL_TET_PAIRS[e].second;
        const size_t gi = tet_indices[li], gj = tet_indices[lj];
        const bool reversed = gi > gj;

        // Always ask in the min -> max direction, so every tet sharing this edge
        // gets the identical answer.
        const size_t lo = reversed ? gj : gi, hi = reversed ? gi : gj;
        const Vector3& p_lo = reversed ? tet_positions[lj] : tet_positions[li];
        const Vector3& p_hi = reversed ? tet_positions[li] : tet_positions[lj];

        if (!reversed) {
            query_edge(lo, hi, p_lo, p_hi, edge_isect_ts[e], edge_isect_normals[e],
                       useRobust, recordNormals);
        } else {
            query_edge(lo, hi, p_lo, p_hi, scratch.ts, scratch.normals, useRobust, recordNormals);
            emit_edge_isect(scratch, /*reversed=*/true, edge_isect_ts[e], edge_isect_normals[e],
                            recordNormals);
        }
    }
}


void InputQueryHandler::query_edge(
    size_t /*global_i*/, size_t /*global_j*/,
    const Vector3& /*pi*/, const Vector3& /*pj*/,
    std::vector<double>& /*out_ts*/,
    std::vector<Vector3>& /*out_normals*/,
    bool /*useRobust*/, bool /*recordNormals*/
){
    throw std::logic_error(
        "this InputQueryHandler does not support single-edge queries; "
        "check supports_edge_query() before using the edge cache.");
}

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

std::vector<std::string> SDFQueryHandler::available_sdf_names() {
    auto names = sdf::getAvailableSDFs();
    std::sort(names.begin(), names.end());
    return names;
}

bool SDFQueryHandler::is_valid_sdf_name(const std::string& name) {
    auto names = sdf::getAvailableSDFs();
    return std::find(names.begin(), names.end(), name) != names.end();
}

SDFQueryHandler::SDFQueryHandler(const std::string& name, float step_size)
    : sdf_name(name), min_step_size(step_size)
{
    // Validate up front so we fail immediately (instead of lazily mid-query).
    if (!is_valid_sdf_name(name))
        throw std::runtime_error("unknown SDF '" + name + "'");

    // Create the SDF evaluation function
    sdf_func = [name](const Vector3& p) -> float {
        glm::vec3 glm_p(p.x, p.y, p.z);
        return sdf::evaluate(name, glm_p);
    };

    log_info("loaded SDF: " + sdf_name);
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
    if (!collect_query_stats) return;
    total_queries += query_count_per_edge[0] + query_count_per_edge[1] + query_count_per_edge[2] + query_count_per_edge[3] + query_count_per_edge[4] + query_count_per_edge[5];
    // Update edge query counts for logging
    local_edge_query_counts = query_count_per_edge;
}


void SDFQueryHandler::query_edge(
    size_t /*global_i*/, size_t /*global_j*/,
    const Vector3& pi, const Vector3& pj,
    std::vector<double>& out_ts,
    std::vector<Vector3>& out_normals,
    bool /*useRobust*/, bool /*recordNormals*/
){
    // Like query_intersections, the SDF march always produces normals.
    std::vector<float> ts_f;
    size_t query_count = 0;
    edge_intersects_SDF(pi, pj, sdf_func, ts_f, out_normals, query_count, min_step_size);

    out_ts.clear();
    out_ts.reserve(ts_f.size());
    for (float t : ts_f) out_ts.push_back(static_cast<double>(t));

    if (collect_query_stats) total_queries += query_count;
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
    log_info("built FCPW accelerator for mesh with " + std::to_string(positions.size()) +
             " vertices and " + std::to_string(polygons.size()) + " faces.");
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

void MeshQueryHandler::query_edge(
    size_t /*global_i*/, size_t /*global_j*/,
    const Vector3& pi, const Vector3& pj,
    std::vector<double>& out_ts,
    std::vector<Vector3>& out_normals,
    bool useRobust, bool recordNormals
){
    find_single_edge_intersections_fcpw(
        pi, pj, accel, out_ts, out_normals, useRobust, recordNormals
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
    double normalize_scale,
    unsigned int seed
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
        std::mt19937 rng(seed);                              // default 1u reproduces the historical offset
        std::normal_distribution<double> gaussian(0.0, 1.0);
        Vector3 dir{gaussian(rng), gaussian(rng), gaussian(rng)};
        double nrm = dir.norm();
        Vector3 t = (nrm > 1e-12 ? dir / nrm : Vector3{1.0, 0.0, 0.0}) * t_mag;
        for (Vector3& p : result.positions) p += t;
    }

    return result;
}
