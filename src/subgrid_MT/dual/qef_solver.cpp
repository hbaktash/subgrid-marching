
#include "subgrid_MT/dual/qef_solver.h"


std::vector<Vector3>
QEF_from_boundary_polygons(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3>& positions,
    const std::vector<Vector3>& normals,
    double regularizer_weight
){
    std::vector<Vector3> dual_points;
    for (const auto& poly: polygons){
        if (poly.size() < 2){
            dual_points.push_back(Vector3::undefined());
            continue;
        }
        std::vector<Eigen::Vector3d> poly_positions, poly_normals;
        Eigen::Vector3d avg_pos = Eigen::Vector3d::Zero();
        for (size_t vidx: poly){
            poly_positions.push_back(vector3_to_eigen(positions[vidx]));
            avg_pos += vector3_to_eigen(positions[vidx]);
            poly_normals.push_back(vector3_to_eigen(normals[vidx]));
        }
        avg_pos /= static_cast<double>(poly.size());
        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        Eigen::Vector3d b = Eigen::Vector3d::Zero();
        for (size_t i = 0; i < poly.size(); i++){
            const Eigen::Vector3d& n = poly_normals[i];
            const Eigen::Vector3d& p = poly_positions[i];
            A += n * n.transpose();
            b += n * n.dot(p);
        }
        A += regularizer_weight * Eigen::Matrix3d::Identity();
        b += regularizer_weight * avg_pos;
        Eigen::Vector3d x = A.ldlt().solve(b);
        dual_points.push_back(eigen_to_vector3(x));
    }
    return dual_points;
}


std::vector<Vector3> get_centroids(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3>& positions
){
    std::vector<Vector3> centroids;
    for (const auto& poly: polygons){
        Vector3 centroid = Vector3::zero();
        for (size_t vidx: poly)
            centroid += positions[vidx];
        centroid /= static_cast<double>(poly.size());
        centroids.push_back(centroid);
    }
    return centroids;
}


Vector3
projected_dual_position(
    const Vector3& dual_p,
    const Vector3& centroid,
    const std::array<Vector3,4>& tet_positions
){
    if (is_inside_tet(dual_p, tet_positions))
        return dual_p;
    double eps = 1e-6;
    Vector3 tet_center = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0;
    Vector3 adjusted_centroid = centroid + 1e-6 * (tet_center - centroid);
    double isect_t;
    segment_tet_intersection(dual_p, adjusted_centroid, tet_positions, isect_t);
    Vector3 dir = adjusted_centroid - dual_p;
    return dual_p + (isect_t + eps) * dir;
}


std::vector<Vector3>
projected_dual_positions(
    const std::vector<Vector3>& dual_positions,
    const std::vector<Vector3>& centroids,
    const std::array<Vector3,4>& tet_positions
){
    std::vector<Vector3> projected_positions;
    for (size_t idx = 0; idx < dual_positions.size(); idx++){
        Vector3 dp = dual_positions[idx];
        if (!dp.isDefined()){
            projected_positions.push_back(dp);
            continue;
        }
        projected_positions.push_back(
            projected_dual_position(dp, centroids[idx], tet_positions)
        );
    }
    return projected_positions;
}


std::vector<Vector3>
dual_positions_from_isect_data(
    const std::array<Vector3,4>& tet_positions,
    const std::vector<std::vector<size_t>>& boundary_polygons,
    const std::vector<Vector3>& isect_positions,
    const std::vector<Vector3>& isect_normals,
    double regularizer_weight,
    bool projected
){
    std::vector<Vector3> dual_positions = QEF_from_boundary_polygons(
        boundary_polygons, isect_positions, isect_normals, regularizer_weight);
    if (projected) {
        std::vector<Vector3> centroids = get_centroids(boundary_polygons, isect_positions);
        dual_positions = projected_dual_positions(dual_positions, centroids, tet_positions);
    }
    return dual_positions;
}
