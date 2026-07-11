#include "subgrid_MT/cv_interpolation.h"
#include "common/triangle_soup.h"


std::vector<Vector3>
comb_curve_to_vec3_curve(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    double bulge
){
    vector<Vector3> poly_positions;
    int n = comb_face.vertices.size();
    for (int idx = 0; idx < n; idx++){
        const auto &cv = comb_face.vertices[idx];
        auto next_cv = comb_face.vertices[(idx + 1) % n];
        Vector3 pos = interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0);
        poly_positions.push_back(pos);
        // add bulge vertex for consecutive vertices on the same edge, except for the last one if it's an open curve
        if (bulge > 0 && next_cv.same_edge_as(cv) && !(comb_face.type == CombFaceType::OPEN && idx == n - 1)){
            Vector3 next_pos = interpolate_comb_vertex(next_cv, tet_positions, edge_isect_ts, -1.0);
            Vector3 interior_face_pos = (2.0 * tet_positions[next_cv.k] + tet_positions[cv.i] + tet_positions[cv.j]) / 4.0;
            Vector3 mid_pos = 0.5 * (pos + next_pos);
            mid_pos = (1 - bulge) * mid_pos + bulge * (interior_face_pos);
            poly_positions.push_back(mid_pos);
        }
    }
    return poly_positions;
}


Vector3
interpolate_comb_vertex(
    const CombinatorialVertex &cv,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts,
    double bulge
){
    int i = cv.i, j = cv.j,
        order = cv.order, n_ij = cv.edge_int;
    int edge_idx = edge_pair_to_index(i, j);
    assert((order <= edge_isect_ts[edge_idx].size()) && (cv.edge_int == edge_isect_ts[edge_idx].size()));

    double t;
    if (i < j){
        t = edge_isect_ts[edge_idx][order - 1]; // convert to 0-based index
    }
    else { // isect_ts are ordered for i < j by convention
        int tmp_order = n_ij - order + 1;
        t = 1 - edge_isect_ts[edge_idx][tmp_order - 1];
    }
    Vector3 pos = (1 - t) * tet_positions[i] + t * tet_positions[j];
    if (!cv.scoop_interior_steiner && !cv.scoop_steiner){
        return pos;
    }
    else {
        if (cv.scoop_interior_steiner && cv.scoop_steiner){
            throw std::invalid_argument("interpolate_comb_vertex: vertex cannot be both scoop_steiner and scoop_interior_steiner");
        }
        double next_t;
        int next_order = order + 1;
        if (i < j){
            next_t = edge_isect_ts[edge_idx][next_order - 1];
        }
        else {
            int tmp_next_order = n_ij - next_order + 1;
            next_t = 1 - edge_isect_ts[edge_idx][tmp_next_order - 1];
        }
        Vector3 next_pos = (1 - next_t) * tet_positions[i] + next_t * tet_positions[j];
        Vector3 mid_pos = 0.5 * (next_pos + pos);
        if (cv.scoop_interior_steiner){
            Vector3 mid_other_edge = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3] - tet_positions[i] - tet_positions[j])/2.0;
            return (1 - bulge) * mid_pos + bulge * mid_other_edge;
        }
        else {
            int k = cv.k;
            Vector3 face_interior_point = (2.0*tet_positions[k] + tet_positions[i] + tet_positions[j]) / 4.0;
            return (1 - bulge) * mid_pos + bulge * face_interior_point;
        }
    }
}


Vector3
interpolate_derived_vertex(
    const DerivedVertex &dv,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts
){
    Vector3 pos = Vector3::zero();
    for (size_t p = 0; p < dv.parents.size(); p++){
        Vector3 parent_pos = interpolate_comb_vertex(dv.parents[p], tet_positions, edge_isect_ts, -1.0);
        pos += dv.weights[p] * parent_pos;
    }
    return pos;
}




// ---- Per-type embed helpers ----

static TriangleSoup
embed_triangle_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts
){
    if (comb_face.vertices.size() != 3)
        throw std::invalid_argument("embed_triangle_face: CombFace of type TRIANGLE must have exactly 3 vertices");
    TriangleSoup soup;
    for (const auto &cv: comb_face.vertices){
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {{0, 1, 2}};
    return soup;
}


static TriangleSoup
embed_quad_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts
){
    if (comb_face.vertices.size() != 4)
        throw std::invalid_argument("embed_quad_face: CombFace of type QUAD02 or QUAD13 must have exactly 4 vertices");
    vector<size_t> tri1 = comb_face.type == CombFaceType::QUAD02 ? vector<size_t>{0, 1, 2} : vector<size_t>{0, 1, 3};
    vector<size_t> tri2 = comb_face.type == CombFaceType::QUAD02 ? vector<size_t>{2, 3, 0} : vector<size_t>{1, 2, 3};
    TriangleSoup soup;
    for (const auto &cv: comb_face.vertices){
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {tri1, tri2};
    return soup;
}


static TriangleSoup
embed_octagon_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts
){
    if (comb_face.vertices.size() != 8)
        throw std::invalid_argument("embed_octagon_face: CombFace of type OCTAGON must have exactly 8 vertices");
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    for (size_t i = 0; i < n; i++){
        soup.vertices.push_back(interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(comb_face.vertices[i]));
        soup.faces.push_back({n, i, (i + 1) % n});
    }
    soup.vertices.push_back(interpolate_derived_vertex(comb_face.derived_vertex, tet_positions, edge_isect_ts));
    if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
    return soup;
}


static TriangleSoup
embed_spiral_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    double bulge
){
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    Vector3 mean_pos = Vector3::zero();
    for (size_t i = 0; i < n; i++){
        Vector3 pos = interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, bulge);
        soup.vertices.push_back(pos);
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(comb_face.vertices[i]));
        mean_pos += pos;
        soup.faces.push_back({n, i, (i + 1) % n});
    }
    mean_pos /= n;
    soup.vertices.push_back(mean_pos);
    if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
    return soup;
}


static TriangleSoup
embed_interior_hexagon_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    double bulge
){
    if (comb_face.vertices.size() != 9)
        throw std::invalid_argument("embed_interior_hexagon_face: CombFace of type INTERIOR_HEXAGON with mid_scoop_vertices=true must have exactly 9 vertices (6 original + 3 mid scoop)");
    TriangleSoup soup;
    for (const auto &cv: comb_face.vertices){
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {
        {8, 2, 5},
        {2, 8, 0}, {2, 0, 1},
        {5, 2, 3}, {5, 3, 4},
        {8, 5, 6}, {8, 6, 7}
    };
    return soup;
}


static Vector3
interior_corner_type_bulge_steiner(
    const CombFace &comb_face,
    const vector<Vector3> &positions,
    const array<Vector3,4> &tet_positions,
    double bulge
){
    CombinatorialVertex cv0 = comb_face.vertices[0],
                        cv2 = comb_face.vertices[2];
    int common_i = (cv0.i == cv2.i || cv0.i == cv2.j) ? cv0.i : (cv0.j == cv2.i || cv0.j == cv2.j) ? cv0.j : -1;
    if (common_i == -1)
        throw std::invalid_argument("embed_corner_type_triangle_face: first two vertices must be on the same edge");
    int l = complement({cv0.i, cv0.j, cv2.i, cv2.j})[0];
    Vector3 mid_pos = 0.5 * (positions[0] + positions[2]);
    Vector3 mid_other_edge = (tet_positions[l] + tet_positions[common_i]) / 2.0;
    return (1 - bulge) * mid_pos + bulge * mid_other_edge;
}


static TriangleSoup
embed_interior_corner_type_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    bool mid_scoop_vertices,
    double bulge
){
    TriangleSoup soup;
    for (const auto &cv: comb_face.vertices){
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    size_t n = soup.vertices.size();
    if (!mid_scoop_vertices){
        for (size_t i = 1; i < n - 1; i++)
            soup.faces.push_back({0, i, i + 1});
    } else {
        Vector3 bulge_steiner = interior_corner_type_bulge_steiner(comb_face, soup.vertices, tet_positions, bulge);
        soup.vertices.push_back(bulge_steiner);
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
        for (size_t i = 0; i < n; i++)
            soup.faces.push_back({n, i, (i + 1) % n});
    }
    return soup;
}


static TriangleSoup
embed_smeared_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    double bulge
){
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    for (const auto &cv: comb_face.vertices){
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    for (size_t i = 2; i < n; i++)
        soup.faces.push_back({1, i, (i + 1) % n});
    return soup;
}


TriangleSoup
embed_comb_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>, 6> &edge_isect_ts,
    bool mid_scoop_vertices,
    double bulge
){
    if (comb_face.vertices.size() < 3)
        throw std::invalid_argument("embed_comb_face: CombFace expects at least 3 vertices. Tunnels should have already received the extra cv's if needed.");

    if (comb_face.type == CombFaceType::TRIANGLE)
        return embed_triangle_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::QUAD02 || comb_face.type == CombFaceType::QUAD13)
        return embed_quad_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::OCTAGON)
        return embed_octagon_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::SPIRAL)
        return embed_spiral_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else if (comb_face.type == CombFaceType::INTERIOR_CORNER_TYPE)
        return embed_interior_corner_type_face(comb_face, tet_positions, edge_isect_ts, mid_scoop_vertices, bulge);
    else if (mid_scoop_vertices && comb_face.type == CombFaceType::INTERIOR_HEXAGON)
        return embed_interior_hexagon_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else if (comb_face.type == CombFaceType::SMEARED_QUAD ||
                comb_face.type == CombFaceType::INTERIOR_HEXAGON ||
                comb_face.type == CombFaceType::TUNNEL_QUAD)
        return embed_smeared_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else
        throw std::invalid_argument("embed_comb_face: CombFace type not supported");
}
