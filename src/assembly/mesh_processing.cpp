#include "assembly/mesh_processing.h"

std::pair<SurfaceMesh*, VertexPositionGeometry*>
process_reconstructed_geometry(
    std::vector<std::vector<size_t>>& face_indices,
    std::vector<Vector3>& vertex_positions,
    bool merge /*= true*/
){
    SimplePolygonMesh ns_simple_mesh(face_indices, vertex_positions);
    if (merge) {
        ns_simple_mesh.mergeIdenticalVertices();
    }
    ns_simple_mesh.stripUnusedVertices();
    ns_simple_mesh.consistentlyOrientFaces();
    SurfaceMesh *ns_mesh = new SurfaceMesh(ns_simple_mesh.polygons);
    ns_mesh->compress();
    ns_mesh->greedilyOrientFaces();
    VertexPositionGeometry* ns_geo = new VertexPositionGeometry(*ns_mesh);
    for (Vertex v : ns_mesh->vertices()){
        ns_geo->inputVertexPositions[v] = ns_simple_mesh.vertexCoordinates[v.getIndex()];
    }
    return {ns_mesh, ns_geo};
}


std::pair<SurfaceMesh*, VertexPositionGeometry*>
mergeIdenticalVertices(
    const double eps,
    std::vector<std::vector<size_t>>& polygons,
    std::vector<Vector3>& vertexCoordinates
) {
    if (eps == 0.0) {
        return process_reconstructed_geometry(polygons, vertexCoordinates);
    }
    const size_t N = vertexCoordinates.size();
    if (N == 0) return {nullptr, nullptr};

    VertexCloudAdaptor cloud{vertexCoordinates};
    KDTree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::vector<size_t> parent(N);
    std::iota(parent.begin(), parent.end(), 0);

    auto find_root = [&](size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](size_t a, size_t b) {
        a = find_root(a);
        b = find_root(b);
        if (a == b) return;
        if (a < b) parent[b] = a;
        else       parent[a] = b;
    };

    const double radius2 = eps * eps;
    using IndexType = KDTree::IndexType;
    using DistanceType = KDTree::DistanceType;
    using ResultItem = nanoflann::ResultItem<IndexType, DistanceType>;

    std::vector<ResultItem> matches;

    for (size_t i = 0; i < N; ++i) {
        double q[3] = { vertexCoordinates[i].x,
                        vertexCoordinates[i].y,
                        vertexCoordinates[i].z };

        matches.clear();
        const size_t nFound = index.radiusSearch(q, radius2, matches);
        (void)nFound;

        for (const auto& m : matches) {
            size_t j = static_cast<size_t>(m.first);
            if (j == i) continue;
            if (j > i) unite(i, j);
        }
    }

    std::vector<size_t> oldToNew(N, std::numeric_limits<size_t>::max());
    std::vector<Vector3> newVerts;
    newVerts.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        size_t r = find_root(i);
        if (oldToNew[r] == std::numeric_limits<size_t>::max()) {
            oldToNew[r] = newVerts.size();
            newVerts.push_back(vertexCoordinates[r]);
        }
        oldToNew[i] = oldToNew[r];
    }

    vertexCoordinates.swap(newVerts);

    for (auto& face : polygons) {
        for (size_t& vi : face) {
            vi = oldToNew[vi];
        }
    }
    return process_reconstructed_geometry(polygons, vertexCoordinates, false);
}


bool check_edge_manifoldness(
    SurfaceMesh& mesh
){
    for (Edge e : mesh.edges()) {
        if (e.halfedge().twin().getIndex() != INVALID_IND &&
            e.halfedge().twin().twin() != e.halfedge()
            ){
            return false;
        }
    }
    return true;
}
