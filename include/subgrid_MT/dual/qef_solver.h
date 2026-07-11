#pragma once

#include "common/utils.h"
#include "query/intersection_query.h"

#include <Eigen/Dense>

inline Eigen::Vector3d vector3_to_eigen(const Vector3& v){
    return Eigen::Vector3d(v.x, v.y, v.z);
}

inline Vector3 eigen_to_vector3(const Eigen::Vector3d& ev){
    return Vector3{ev.x(), ev.y(), ev.z()};
}


// NOTE: convert_mod2_normals moved to query/mod2_reduction.h

std::vector<Vector3>
QEF_from_boundary_polygons(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3>& positions,
    const std::vector<Vector3>& normals,
    double regularizer_weight = 1e-3
);

std::vector<Vector3> get_centroids(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3>& positions
);

Vector3
projected_dual_position(
    const Vector3& pd,
    const Vector3& centroid,
    const std::array<Vector3,4>& tet_positions
);

std::vector<Vector3>
projected_dual_positions(
    const std::vector<Vector3>& dual_positions,
    const std::vector<Vector3>& centroids,
    const std::array<Vector3,4>& tet_positions
);

std::vector<Vector3>
dual_positions_from_isect_data(
    const std::array<Vector3,4>& tet_positions,
    const std::vector<std::vector<size_t>>& boundary_polygons,
    const std::vector<Vector3>& isect_positions,
    const std::vector<Vector3>& isect_normals,
    double regularizer_weight,
    bool projected
);
