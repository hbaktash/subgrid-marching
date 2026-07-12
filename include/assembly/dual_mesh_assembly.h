#pragma once

#include "assembly/mesh_processing.h"

void cleanup_dual_polygons(
    const std::vector<std::vector<size_t>>& input_dual_polygons,
    const std::vector<Vector3>& input_dual_positions,
    std::vector<std::vector<size_t>>& output_dual_polygons,
    std::vector<Vector3>& output_dual_positions
);

std::vector<size_t>
dual_polygon_manifold_unori(
    Vertex v
);

std::vector<std::vector<Face>>
neigh_face_components(
    Vertex v
);

void
build_primal_twin_halfedge_map_from_face_per_edge_data(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<std::vector<std::array<size_t, 3>>>& tet_face_indices,
    std::vector<std::vector<std::pair<size_t, size_t>>>& twins
);

void remove_pinch_vertices(
    SurfaceMesh& mesh,
    VertexPositionGeometry& geo,
    const std::vector<std::vector<size_t>>& in_polygons,
    std::vector<std::vector<size_t>>& out_polygons,
    std::vector<Vector3>& out_positions
);

std::pair<SurfaceMesh*, VertexPositionGeometry*>
construct_primal_mesh_from_face_per_edge_data(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3> &positions,
    const std::vector<std::vector<std::array<size_t, 3>>>& tet_face_indices
);

std::pair<SurfaceMesh*, VertexPositionGeometry*>
construct_dual_mesh(
    SurfaceMesh& primal_mesh,
    VertexPositionGeometry& primal_geo,
    const std::vector<Vector3>& dual_positions
);

