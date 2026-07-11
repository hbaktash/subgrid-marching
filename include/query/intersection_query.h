#pragma once

#include "common/utils.h"
#include <memory>

#include "fcpw/fcpw.h"
#include "query/fcpw_robust_queries.h"


// NOTE: convert_mod2_intersections moved to query/mod2_reduction.h

inline std::array<int,6> counts_from_isect_ts(const std::array<std::vector<double>,6>& isect_ts) {
    std::array<int,6> counts;
    for (size_t e = 0; e < 6; ++e) counts[e] = (int)isect_ts[e].size();
    return counts;
}

// single edge-tet intersection query
bool segment_tet_intersection(
    const geometrycentral::Vector3& p0,
    const geometrycentral::Vector3& p1,
    const std::array<geometrycentral::Vector3,4>& tet_positions,
    double& out_t
); 


std::array<std::vector<double>, 6>
two_fold_single_tet_edge_intersections(
    const std::array<std::vector<double>, 6>& base_intersections
);

std::array<std::vector<geometrycentral::Vector3>, 6>
two_fold_single_tet_edge_isect_normals(
    const std::array<std::vector<geometrycentral::Vector3>, 6>& isect_normals
);

// Build an accelerated intersection structure for the input surface mesh using FCPW.
// The resulting fcpw::Scene is consumed by find_single_tet_edge_intersections_fcpw
// and the closest-point queries below (HAVE_FCPW must be defined at compile-time).
void
build_fcpw_accelerator(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<geometrycentral::Vector3>& vertexCoordinates,
    fcpw::Scene<3>& out_accel
);

void build_fcpw_accelerator(
    SurfaceMesh& input_mesh,
    VertexPositionGeometry& input_geometry,
    fcpw::Scene<3>& out_accel
);


// Build FCPW accelerator from polygon soup
fcpw::Scene<3>
get_fcpw_accel(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<geometrycentral::Vector3>& vertexCoordinates
);

fcpw::Scene<3>
get_fcpw_accel(
    SurfaceMesh& input_mesh,
    VertexPositionGeometry& input_geometry
);

// use FCPW for
//   - tetMeshEdges-mesh intersection query
//   -closest point queries

void
find_single_tet_edge_intersections_fcpw(
    const std::array<geometrycentral::Vector3,4>& tet_positions,
    fcpw::Scene<3>& accel,
    std::array<std::vector<double>, 6>& out_edge_isect_ts,
    std::array<std::vector<geometrycentral::Vector3>, 6>& out_isect_normals,
    bool useRobust = false,
    bool recordNormals = true
);



// -------------------------------- CP queries; with FCPW --------------------------------

void closest_point_fcpw(
    const Vector3 &query_point,
    fcpw::Scene<3>& accel,
    geometrycentral::Vector3 &out_closest_point,
    geometrycentral::Vector3 &out_closest_normal
);