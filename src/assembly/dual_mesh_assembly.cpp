#include "assembly/dual_mesh_assembly.h"

void cleanup_dual_polygons(
    const std::vector<std::vector<size_t>>& input_dual_polygons,
    const std::vector<Vector3>& input_dual_positions,
    std::vector<std::vector<size_t>>& output_dual_polygons,
    std::vector<Vector3>& output_dual_positions
){
    output_dual_polygons.clear();
    std::vector<bool> vertex_used(input_dual_positions.size(), false);
    size_t removed_count = 0;
    for (const auto &poly : input_dual_polygons){
        if (poly.size() < 3){
            removed_count++;
            continue;
        }
        output_dual_polygons.push_back(poly);
        for (auto vidx : poly){
            vertex_used[vidx] = true;
        }
    }
    if (removed_count > 0)
        log_info("cleanup_dual_polygons: removed " + std::to_string(removed_count) + "/" +
                 std::to_string(input_dual_polygons.size()) + " degenerate (<3) polygons.");
    std::vector<size_t> old_to_new_index(input_dual_positions.size(), SIZE_MAX);
    output_dual_positions.clear();
    for (size_t vidx = 0; vidx < input_dual_positions.size(); ++vidx){
        if (vertex_used[vidx]){
            size_t new_idx = output_dual_positions.size();
            old_to_new_index[vidx] = new_idx;
            output_dual_positions.push_back(input_dual_positions[vidx]);
        }
    }
    for (auto &poly : output_dual_polygons){
        for (auto &vidx : poly){
            vidx = old_to_new_index[vidx];
        }
    }
}


std::vector<size_t>
dual_polygon_manifold_unori(
    Vertex v
){
    Edge start_e = v.halfedge().edge();
    for (Edge e: v.adjacentEdges()){
        if (e.isBoundary()){
            start_e = e;
            if (e.halfedge().twin() != e.halfedge()){
                throw std::runtime_error("dual_polygon_manifold_unori: boundary edge with non-boundary twin detected.");
            }
            break;
        }
    }

    std::vector<size_t> dualpoly;
    Halfedge curr_he = start_e.halfedge();
    Halfedge start_he = curr_he;
    do {
        Face f = curr_he.face();
        if (f.isBoundaryLoop()){
            throw std::runtime_error("dual_polygon_manifold_unori: boundary loop face detected.");
        }
        dualpoly.push_back(f.getIndex());
        Halfedge next_he;
        if (curr_he.tipVertex() == v){
            next_he = curr_he.next().twin();
        } else if (curr_he.tailVertex() == v){
            Halfedge f_prev_he = curr_he.prevOrbitFace();
            if (f_prev_he.next() != curr_he){
                throw std::runtime_error("dual_polygon_manifold_unori: previous HE didn't match.");
            }
            next_he = f_prev_he.twin();
        } else {
            throw std::runtime_error("dual_polygon_manifold_unori: halfedge not incident to vertex.");
        }
        curr_he = next_he;
    } while (curr_he != start_he);
    return dualpoly;
}

std::vector<std::vector<Face>>
neigh_face_components(
    Vertex v
){
    std::vector<std::vector<Face>> components;
    std::unordered_set<Face> seen;

    for (Face seed_face : v.adjacentFaces()) {
        if (seen.find(seed_face) != seen.end()) continue;

        std::vector<Face> component;
        std::vector<Face> to_proc{seed_face};

        while (!to_proc.empty()) {
            Face f = to_proc.back();
            to_proc.pop_back();

            if (seen.find(f) != seen.end()) continue;
            seen.insert(f);
            component.push_back(f);

            for (Edge e : f.adjacentEdges()) {
                if (e.firstVertex() != v && e.secondVertex() != v) continue;
                for (Face fn : e.adjacentFaces()) {
                    if (seen.find(fn) == seen.end() && fn != f) {
                        to_proc.push_back(fn);
                    }
                }
            }
        }
        if (!component.empty()) {
            components.push_back(component);
        }
    }
    return components;
}


void
build_primal_twin_halfedge_map_from_face_per_edge_data(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<std::vector<std::array<size_t, 3>>>& tet_face_indices,
    std::vector<std::vector<std::pair<size_t, size_t>>>& twins
){
    if (polygons.size() != tet_face_indices.size()){
        throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: polygons and tet_face_indices size mismatch.");
    }

    size_t nFaces = polygons.size();
    twins.clear();
    twins.resize(nFaces);

    for (size_t f = 0; f < nFaces; ++f) {
        size_t deg = polygons[f].size();
        if (tet_face_indices[f].size() != deg) {
            throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: per-face tet_face_indices size does not match polygon degree.");
        }
        twins[f].assign(deg, {INVALID_IND, INVALID_IND});
    }

    struct HalfedgeRef {
        size_t faceIdx;
        size_t edgeIdx;
    };

    struct EdgeFaceKey {
        size_t a, b;
        std::array<size_t, 3> face;

        bool operator<(const EdgeFaceKey& other) const {
            if (a < other.a) return true;
            if (a > other.a) return false;
            if (b < other.b) return true;
            if (b > other.b) return false;
            for (int i = 0; i < 3; ++i) {
                if (face[i] < other.face[i]) return true;
                if (face[i] > other.face[i]) return false;
            }
            return false;
        }
    };

    std::map<EdgeFaceKey, std::vector<HalfedgeRef>> buckets;

    for (size_t f = 0; f < nFaces; ++f) {
        const auto& poly = polygons[f];
        const auto& faceTriplets = tet_face_indices[f];
        size_t deg = poly.size();
        if (deg < 2) {
            throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: degenerate polygon with less than 2 vertices.");
        }

        for (size_t e = 0; e < deg; ++e) {
            size_t v0 = poly[e];
            size_t v1 = poly[(e + 1) % deg];
            if (v0 == v1) {
                throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: zero-length edge found in polygon.");
            }

            EdgeFaceKey key;
            key.a = std::min(v0, v1);
            key.b = std::max(v0, v1);

            key.face = faceTriplets[e];
            std::sort(key.face.begin(), key.face.end());

            buckets[key].push_back({f, e});
        }
    }

    for (auto& kv : buckets) {
        auto& refs = kv.second;
        if (refs.size() <= 1) continue;
        if (refs.size() > 2) {
            throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: more than two halfedges share the same edge and face.");
        }

        const auto& h0 = refs[0];
        const auto& h1 = refs[1];

        if (h0.faceIdx == h1.faceIdx) {
            throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: two halfedges of the same polygon share the same edge and face.");
        }

        if (twins[h0.faceIdx][h0.edgeIdx].first != INVALID_IND ||
            twins[h1.faceIdx][h1.edgeIdx].first != INVALID_IND) {
            throw std::runtime_error("build_primal_twin_halfedge_map_from_face_per_edge_data: halfedge already has a twin when attempting to assign a new one.");
        }

        twins[h0.faceIdx][h0.edgeIdx] = {h1.faceIdx, h1.edgeIdx};
        twins[h1.faceIdx][h1.edgeIdx] = {h0.faceIdx, h0.edgeIdx};
    }

    // Remaining INVALID_IND entries are boundary: self-twin for SurfaceMesh constructor
    for (size_t f = 0; f < nFaces; ++f) {
        size_t deg = polygons[f].size();
        for (size_t e = 0; e < deg; ++e) {
            if (twins[f][e].first == INVALID_IND) {
                twins[f][e] = {f, e};
            }
        }
    }
}


void remove_pinch_vertices(
    SurfaceMesh& mesh,
    VertexPositionGeometry& geo,
    const std::vector<std::vector<size_t>>& in_polygons,
    std::vector<std::vector<size_t>>& out_polygons,
    std::vector<Vector3>& out_positions
){
    size_t nVerts = mesh.nVertices();
    out_positions = std::vector<Vector3>(nVerts);
    for (Vertex v : mesh.vertices()){
        out_positions[v.getIndex()] = geo.inputVertexPositions[v];
    }

    FaceData<std::vector<std::pair<size_t, size_t>>> pinch_edits(mesh);
    for (Vertex v : mesh.vertices()){
        std::vector<std::vector<Face>> neigh_face_comps = neigh_face_components(v);
        if (neigh_face_comps.size() > 1){
            for (size_t i = 1; i < neigh_face_comps.size(); ++i){
                std::vector<Face> comp = neigh_face_comps[i];
                size_t new_v_id = nVerts;
                out_positions.push_back(geo.inputVertexPositions[v]);
                nVerts++;
                std::pair<size_t, size_t> idx_change = {v.getIndex(), new_v_id};
                for (Face f : comp){
                    pinch_edits[f].push_back(idx_change);
                }
            }
        }
    }

    out_polygons = std::vector<std::vector<size_t>>(in_polygons.size());
    for (Face f: mesh.faces()){
        std::vector<size_t> new_poly = in_polygons[f.getIndex()];
        for (auto idx_change : pinch_edits[f]){
            size_t old_idx = idx_change.first;
            size_t new_idx = idx_change.second;
            for (size_t &vidx : new_poly){
                if (vidx == old_idx){
                    vidx = new_idx;
                }
            }
        }
        out_polygons[f.getIndex()] = new_poly;
    }
}


std::pair<SurfaceMesh*, VertexPositionGeometry*>
construct_primal_mesh_from_face_per_edge_data(
    const std::vector<std::vector<size_t>>& polygons,
    const std::vector<Vector3> &positions,
    const std::vector<std::vector<std::array<size_t, 3>>>& tet_face_indices
){
    std::vector<std::vector<std::pair<size_t, size_t>>> twins;
    build_primal_twin_halfedge_map_from_face_per_edge_data(polygons, tet_face_indices, twins);
    std::vector<std::vector<std::tuple<size_t, size_t>>> twins_tuple;
    pairss_to_tupless(twins, twins_tuple);
    SurfaceMesh* initial_primal_mesh = new SurfaceMesh(polygons, twins_tuple);
    VertexPositionGeometry *initial_geometry = new VertexPositionGeometry(*initial_primal_mesh);
    for (Vertex v : initial_primal_mesh->vertices()){
        initial_geometry->inputVertexPositions[v] = positions[v.getIndex()];
    }

    if (!check_edge_manifoldness(*initial_primal_mesh)){
        throw std::runtime_error("construct_primal_mesh_from_face_per_edge_data: mesh is not edge-manifold.");
    }

    bool has_pinch_vertex = false;
    for (Vertex v : initial_primal_mesh->vertices()){
        if (neigh_face_components(v).size() > 1){
            has_pinch_vertex = true;
            break;
        }
    }
    if (!has_pinch_vertex){
        initial_primal_mesh->greedilyOrientFaces();
        return {initial_primal_mesh, initial_geometry};
    }

    log_warn("pinch vertices detected in the primal mesh; separating the pinches.");
    std::vector<std::vector<size_t>> out_polygons;
    std::vector<Vector3> out_positions;
    remove_pinch_vertices(*initial_primal_mesh, *initial_geometry, polygons, out_polygons, out_positions);
    delete initial_geometry;
    delete initial_primal_mesh;

    SurfaceMesh *fixed_mesh = new SurfaceMesh(out_polygons, twins_tuple);
    fixed_mesh->greedilyOrientFaces();
    VertexPositionGeometry* fixed_geometry = new VertexPositionGeometry(*fixed_mesh);
    for (Vertex v : fixed_mesh->vertices()){
        fixed_geometry->inputVertexPositions[v] = out_positions[v.getIndex()];
    }
    return {fixed_mesh, fixed_geometry};
}

std::pair<SurfaceMesh*, VertexPositionGeometry*>
construct_dual_mesh(
    SurfaceMesh& primal_mesh,
    VertexPositionGeometry& primal_geo,
    const std::vector<Vector3>& dual_positions
){
    std::vector<std::vector<size_t>> dual_polygons(primal_mesh.nVertices());
    for (Vertex v : primal_mesh.vertices()){
        dual_polygons[v.getIndex()] = dual_polygon_manifold_unori(v);
    }

    std::vector<std::vector<size_t>> triangulated_dual_polygons;
    triangulate_polygons(dual_polygons, triangulated_dual_polygons, true);

    std::vector<std::vector<size_t>> cleaned_dual_polygons;
    std::vector<Vector3> cleaned_dual_positions;
    cleanup_dual_polygons(triangulated_dual_polygons, dual_positions,
                    cleaned_dual_polygons, cleaned_dual_positions);

    SurfaceMesh* dual_mesh = new SurfaceMesh(cleaned_dual_polygons);
    VertexPositionGeometry* dual_geo = new VertexPositionGeometry(*dual_mesh);
    for (Vertex dv : dual_mesh->vertices()){
        dual_geo->inputVertexPositions[dv] = cleaned_dual_positions[dv.getIndex()];
    }
    dual_mesh->greedilyOrientFaces();
    return {dual_mesh, dual_geo};
}
