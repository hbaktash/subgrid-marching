#include "query/intersection_query.h"


// single edge-tet intersection query
bool segment_tet_intersection(
    const geometrycentral::Vector3& p0,
    const geometrycentral::Vector3& p1,
    const std::array<geometrycentral::Vector3,4>& tet_positions,
    double& out_t
){
    // intersect for each face and return the closest intersection
    out_t = std::numeric_limits<double>::infinity();
    for (auto triplet: ALL_TET_TRIPLETS){
        Vector3 a = tet_positions[triplet[0]];
        Vector3 b = tet_positions[triplet[1]];
        Vector3 c = tet_positions[triplet[2]];
        double t_face;
        if (segment_triangle_intersection(p0, p1, a, b, c, t_face)){
            out_t = std::min(out_t, t_face);
        }
    }
    return out_t != std::numeric_limits<double>::infinity();
} 



std::array<std::vector<double>, 6>
two_fold_single_tet_edge_intersections(
    const std::array<std::vector<double>, 6>& base_intersections
){
    std::array<std::vector<double>, 6> two_fold_intersections;
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        std::vector<double> two_fold_ts;
        for (double t : base_intersections[eIdx]) {
            two_fold_ts.push_back(t);
            two_fold_ts.push_back(t); // duplicate
        }
        two_fold_intersections[eIdx] = two_fold_ts;
    }
    return two_fold_intersections;
}


std::array<std::vector<geometrycentral::Vector3>, 6>
two_fold_single_tet_edge_isect_normals(
    const std::array<std::vector<geometrycentral::Vector3>, 6>& isect_normals
){
    std::array<std::vector<geometrycentral::Vector3>, 6> two_fold_normals;
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        std::vector<geometrycentral::Vector3> two_fold_ns;
        for (const geometrycentral::Vector3& n : isect_normals[eIdx]) {
            two_fold_ns.push_back(n);
            two_fold_ns.push_back(n); // duplicate
        }
        two_fold_normals[eIdx] = two_fold_ns;
    }
    return two_fold_normals;
}


#ifdef HAVE_FCPW

// Build an FCPW mesh + BVH accelerator. We store the BVH (or mesh+bvh) as an opaque pointer in out_accel.

void
build_fcpw_accelerator(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<geometrycentral::Vector3>& vertexCoordinates,
    fcpw::Scene<3>& out_accel
){
    // turn inout positions and faces into nx3 matrix
    Eigen::MatrixXf pos_mat(vertexCoordinates.size(), 3);
    for (size_t i = 0; i < vertexCoordinates.size(); ++i) {
        Vector3 pos = vertexCoordinates[i];
        pos_mat(i, 0) = pos.x;
        pos_mat(i, 1) = pos.y;
        pos_mat(i, 2) = pos.z;
    }

    Eigen::MatrixXi face_mat(polygons.size(), 3);
    for (size_t i = 0; i < polygons.size(); ++i) {
        auto face = polygons[i];
        if (face.size() != 3) 
            throw std::logic_error("Only triangular faces are supported; call triangulate");
        face_mat(i, 0) = face[0];
        face_mat(i, 1) = face[1];
        face_mat(i, 2) = face[2];
    }

    // initialize a 3d scene
    fcpw::Scene<3> scene;

    // load positions and indices of a single triangle mesh
    scene.setObjectCount(1);
    scene.setObjectVertices(pos_mat, 0);
    scene.setObjectTriangles(face_mat, 0);

    // build acceleration structure
    fcpw::AggregateType aggregateType = fcpw::AggregateType::Bvh_SurfaceArea;
    bool buildVectorizedBvh = true;
    scene.build(aggregateType, buildVectorizedBvh);

    out_accel = std::move(scene);
    // out_accel = scene;
}

void build_fcpw_accelerator(
    SurfaceMesh& input_mesh,
    VertexPositionGeometry& input_geometry,
    fcpw::Scene<3>& out_accel
){
    std::vector<std::vector<size_t>> polygons = input_mesh.getFaceVertexList();
    std::vector<geometrycentral::Vector3> vertexCoordinates(input_mesh.nVertices());
    for (Vertex v : input_mesh.vertices()) {
        vertexCoordinates[v.getIndex()] = input_geometry.vertexPositions[v];
    }
    build_fcpw_accelerator(polygons, vertexCoordinates, out_accel);
}

fcpw::Scene<3>
get_fcpw_accel(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<geometrycentral::Vector3>& vertexCoordinates
){
    // build FCPW accelerator
    fcpw::Scene<3> local_accel;
    build_fcpw_accelerator(polygons, vertexCoordinates, local_accel);
    return local_accel;
}

fcpw::Scene<3>
get_fcpw_accel(
    SurfaceMesh& input_mesh,
    VertexPositionGeometry& input_geometry
){
    // build FCPW accelerator
    fcpw::Scene<3> local_accel;
    build_fcpw_accelerator(input_mesh, input_geometry, local_accel);
    return local_accel;
}

void
find_single_tet_edge_intersections_fcpw(
    const std::array<geometrycentral::Vector3,4>& tet_positions,
    fcpw::Scene<3>& accel,
    std::array<std::vector<double>, 6>& out_edge_isect_ts,
    std::array<std::vector<geometrycentral::Vector3>, 6>& out_isect_normals,
    bool useRobust,
    bool recordNormals
){
    // intersect with all edges of vol_mesh
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        std::pair<int,int> edge = ALL_TET_PAIRS[eIdx];
        Vector3 a = tet_positions[edge.first],
                b = tet_positions[edge.second];
        fcpw::Vector3 origin;
        origin << a.x, a.y, a.z;
        fcpw::Vector3 dir;
        dir << b.x - a.x, b.y - a.y, b.z - a.z;
        double seg_len = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        // if (seg_len <= 1e-8) continue;
        for (int k=0;k<3;++k) dir[k] /= seg_len;
        // intersect using FCPW
        fcpw::Ray<3> ray(origin, dir, seg_len + (useRobust ? 0. : 1e-6)); // max length 1.0
        std::vector<fcpw::Interaction<3>> isects;
        if (useRobust) {
            isects = intersect_robust_all_hits(accel, ray, recordNormals, 1e-6, -1);
        } else {
            accel.intersect(ray, isects, false, true);
        }

        // convert intersections to t parameters along segment
        for (size_t ii = 0; ii < isects.size(); ii++) {
            const fcpw::Interaction<3>& interaction = isects[ii];
            if (interaction.d < 0 || interaction.d > seg_len) {
                continue;
            }
            double t_param = interaction.d / seg_len;
            out_edge_isect_ts[eIdx].push_back(t_param);
            if (recordNormals) {
                out_isect_normals[eIdx].push_back(geometrycentral::Vector3{interaction.n[0], interaction.n[1], interaction.n[2]});
            }
        }
    }
}


// --------------------  Closest point queries using FCPW --------------------

void closest_point_fcpw(
    const Vector3 &query_point,
    fcpw::Scene<3>& accel,
    geometrycentral::Vector3 &out_closest_point,
    geometrycentral::Vector3 &out_closest_normal
){
    fcpw::Interaction<3> interaction;
    fcpw::Vector3 query;
    query << query_point.x, query_point.y, query_point.z;
    accel.findClosestPoint(query, interaction, 1e10, true); // 
    fcpw::Vector3 p = interaction.p;
    out_closest_point = geometrycentral::Vector3{p[0], p[1], p[2]};
    fcpw::Vector3 n = interaction.n;
    out_closest_normal = geometrycentral::Vector3{n[0], n[1], n[2]};
}


// NOTE: the !HAVE_FCPW fallback block was removed during clean-up. It was
// already non-functional (input_query_handler unconditionally calls the
// FCPW-only entry points).

#endif