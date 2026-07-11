#pragma once

#include "common/utils.h"

#include "nanoflann.hpp"

std::pair<SurfaceMesh*, VertexPositionGeometry*>
process_reconstructed_geometry(
    std::vector<std::vector<size_t>>& face_indices,
    std::vector<Vector3>& vertex_positions,
    bool merge = true
);


// KDTree adaptor for epsilon-based vertex merging
struct VertexCloudAdaptor {
    const std::vector<Vector3>& pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline double kdtree_distance(const double* p, const size_t idx, size_t /*dim*/) const {
        const Vector3& q = pts[idx];
        double dx = p[0] - q.x;
        double dy = p[1] - q.y;
        double dz = p[2] - q.z;
        return dx*dx + dy*dy + dz*dz;
    }

    inline double kdtree_get_pt(const size_t idx, int dim) const {
        if (dim == 0) return pts[idx].x;
        if (dim == 1) return pts[idx].y;
        return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, VertexCloudAdaptor>,
    VertexCloudAdaptor,
    3
>;

bool check_edge_manifoldness(
    SurfaceMesh& mesh
);

std::pair<SurfaceMesh*, VertexPositionGeometry*>
mergeIdenticalVertices(
    const double eps,
    std::vector<std::vector<size_t>>& polygons,
    std::vector<Vector3>& vertexCoordinates
);
