#include "subgrid_MT/evensum_construction.h"
#include "subgrid_MT/simplicial_embedding.h"



static int
get_interior_vertex(
    const vector<int> &separated_vs
){
    int interior_vertex;
    if (separated_vs.size() == 4 || separated_vs.size() == 0){
        interior_vertex = -1;
    }
    else if (separated_vs.size() == 3){
        vector<int> complement_vs = complement(separated_vs);
        interior_vertex = complement_vs[0];
    }
    else if (separated_vs.size() == 1){
        interior_vertex = separated_vs[0];
    }
    else {
        throw std::logic_error("get_interior_vertex: Invalid number of separated vertices");
    }
    if (interior_vertex < -1 || interior_vertex >= 4){
        throw std::logic_error("get_interior_vertex: computed interior vertex out of range [-1, 3]");
    }
    return interior_vertex;
}

EdgeInts
combface_to_edgeints(
    const CombFace &comb_face
){
    EdgeInts result;
    for (auto &cv: comb_face.vertices){
        int i = cv.i, j = cv.j;
        result[{i,j}]++;
    }
    return result;
}


vector<int>
separated_vertices(
    const EdgeInts &cruve_edge_ints
){
    // include vertex 0 by default
    // any other vertex is added depending on edge (0,i) count parity
    vector<int> result = {0};
    for (int i = 1; i < 4; i++){
        if (cruve_edge_ints[{0,i}] % 2 == 0){
            result.push_back(i);
        }
    }
    return result;
}


static pair<array<vector<double>,6>, EdgeOccupations>
extract_curve_edge_isect_ts(
    const vector<CombinatorialVertex> &curve_cvs,
    const array<vector<double>,6> &global_edge_isect_ts
){
    EdgeInts curve_edge_ints(global_edge_isect_ts);
    EdgeOccupations curve_edge_occupations(curve_edge_ints);
    for (const auto &cv: curve_cvs){
        curve_edge_occupations.set(cv, true);
    }
    array<vector<double>,6> curve_edge_isect_ts;
    for (int e = 0; e < 6; e++){
        int i = ALL_TET_PAIRS[e].first, j = ALL_TET_PAIRS[e].second;
        for (size_t k = 0; k < global_edge_isect_ts[e].size(); k++){
            if (curve_edge_occupations.get(i,j,k)){
                curve_edge_isect_ts[e].push_back(global_edge_isect_ts[e][k]);
            }
        }
    }
    return {curve_edge_isect_ts, curve_edge_occupations};
}


static TriangleSoup
emit_comb_face(
    const CombFace &comb_face,
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &edge_isect_ts,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
){
    return embed_comb_face(comb_face, tet_positions, edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge);
}


TriangleSoup
evensum_surface_construction(
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &global_edge_isect_ts,
    const vector<CombFace> &non_normal_comb_curves,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
){
    TriangleSoup soup;
    for (auto &non_normal_comb_curve: non_normal_comb_curves){
        EdgeInts evensum_edge_ints = combface_to_edgeints(non_normal_comb_curve);
        vector<int> separated_vs = separated_vertices(evensum_edge_ints);
        if (separated_vs.size() == 2){ // diagonal type
            soup.add_local_soup(diagonal_surface_construction(tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge));
        }
        else { // non-diagonal curve; corner or contractible type
            soup.add_local_soup(nondiagonal_surface_construction(tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge));
        }
    }
    return soup;
}


// DIAGONAL curve: emit it directly as a spiral face, using GLOBAL indices and isects.
TriangleSoup
diagonal_surface_construction(
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &global_edge_isect_ts,
    const CombFace &diagonal_comb_curve,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
){
    CombFace spiral_comb_face = diagonal_comb_curve;
    spiral_comb_face.type = CombFaceType::SPIRAL;
    if (scoop_mid_vertices){
        spiral_comb_face = add_mid_scoop_vertices_spiral(spiral_comb_face);
    }
    return emit_comb_face(spiral_comb_face, tet_positions, global_edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge);
}


TriangleSoup
nondiagonal_surface_construction(
    const array<Vector3,4> &tet_positions,
    const array<vector<double>,6> &global_edge_isect_ts,
    const CombFace &nondiagonal_comb_curve,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
){
    auto [local_edge_isect_ts, curve_edge_occupations] = extract_curve_edge_isect_ts(nondiagonal_comb_curve.vertices, global_edge_isect_ts);
    EdgeInts smeared_curve_edge_ints = combface_to_edgeints(nondiagonal_comb_curve);
    int interior_vertex = get_interior_vertex(separated_vertices(smeared_curve_edge_ints));
    assert(interior_vertex >= 0 || (smeared_curve_edge_ints[{0,1}] % 2 == 0 && smeared_curve_edge_ints[{0,2}] % 2 == 0 && smeared_curve_edge_ints[{0,3}] % 2 == 0 && smeared_curve_edge_ints[{1,2}] % 2 == 0 && smeared_curve_edge_ints[{1,3}] % 2 == 0 && smeared_curve_edge_ints[{2,3}] % 2 == 0));

    NormalCoordinates smeared_curve_nc(smeared_curve_edge_ints);
    if (smeared_curve_nc.triangle_ineq || !smeared_curve_nc.even_sum) {
        throw std::invalid_argument("smeared_surface_construction: should satisfy even sum and NOT satisfy triangle inequality");
    }

    vector<CombFace> comb_faces;
    auto append = [&comb_faces](const vector<CombFace> &section){ comb_faces.insert(comb_faces.end(), section.begin(), section.end()); };
    append(smeared_quads_construction(smeared_curve_nc, interior_vertex));
    append(interior_region_construction(smeared_curve_nc, interior_vertex));
    if (scoop_mid_vertices){
        append(tunnel_construction(smeared_curve_nc, interior_vertex));
    }

    TriangleSoup curve_soup;
    for (CombFace face : comb_faces){
        if (scoop_mid_vertices){
            face = add_mid_scoop_vertices(face, smeared_curve_nc);
        }
        if (face.vertices.size() < 3){
            continue;
        }
        curve_soup.add_local_soup(emit_comb_face(face, tet_positions, local_edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge));
    }
    if (TriangleSoup::COMB_MERGE)
        curve_soup.remap_orders(global_edge_isect_ts, curve_edge_occupations, true);
    return curve_soup;
}



// SMEARED QUADs across all faces: for each face corner, build the polygons
// behind the even/odd smeared segments (segment_order 1 is the corner triangle).
vector<CombFace>
smeared_quads_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex
){
    vector<CombFace> faces;
    for (const auto &ijk: ALL_TET_TRIPLETS){
        // first go over the corner polygons/quad ,
        // then the interior region polygon
        for (int i: ijk){ // at corner i
            int j = ijk[0] == i ? ijk[1] : ijk[0];
            int k = ijk[0] + ijk[1] + ijk[2] - i - j;
            int n_ij = smeared_curve_nc.edge_ints[{i,j}],
                n_ik = smeared_curve_nc.edge_ints[{i,k}];

            bool i_is_interior = (interior_vertex == i);
            for (int segment_order = 1; segment_order <= smeared_curve_nc.scc.get(i, j, k); segment_order++){ // order is 1-based
                // always considering the polygon behind the current segment
                bool construct = (i_is_interior && segment_order % 2 == 1) || (!i_is_interior && segment_order % 2 == 0);
                if (construct && segment_order > 1){ // segment_order 1 is just the corner triangle, which is handled separately below
                    // construct the polygon for this segment; the cv.k's here are helpers for mid scoop vertex later
                    CombinatorialVertex cvij_1(i, j, segment_order-1, n_ij, k);
                    CombinatorialVertex cvij_2(i, j, segment_order, n_ij, k);
                    CombinatorialVertex cvik_2(i, k, segment_order, n_ik, j);
                    CombinatorialVertex cvik_1(i, k, segment_order-1, n_ik, j);
                    // order matters here for future scoop mid vertex construction
                    faces.push_back(CombFace({cvij_1, cvij_2, cvik_2, cvik_1}, CombFaceType::SMEARED_QUAD));
                }
            }
        }
    }
    return faces;
}


// TUNNELS connecting the two sides of each interior segment (scoop mode only).
// Returns the raw two-vertex tunnel quads; the mid scoop vertices that expand
// them (and the drop of empty ones) are handled by the caller.
vector<CombFace>
tunnel_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex
){
    vector<CombFace> faces;
    for (const auto &ij: ALL_TET_PAIRS){
        int i = ij.first, j = ij.second;
        int n_ij = smeared_curve_nc.edge_ints[{i,j}];
        for (int segment_order = 2; segment_order <= n_ij; segment_order++){ // order is 1-based
            // consider the segment between segment_order-1 and segment_order
            bool i_is_interior = (interior_vertex == i);
            bool construct = (i_is_interior && segment_order % 2 == 1) || (!i_is_interior && segment_order % 2 == 0);
            // build the tunnel quad for this segment
            if (construct){
                // first order vertices
                CombinatorialVertex cvij_1(i, j, segment_order-1, n_ij);
                CombinatorialVertex cvij_2(i, j, segment_order, n_ij);
                faces.push_back(CombFace({cvij_1, cvij_2}, CombFaceType::TUNNEL_QUAD));
            }
        }
    }
    return faces;
}


vector<CombFace>
interior_region_construction(
    const NormalCoordinates &smeared_curve_nc,
    const int &interior_vertex // vertex that is inside the curve; -1 if contractible curve
){
    if (interior_vertex == -1){
        // contractible curve; find the interior region on some face
        for (auto &ijk: ALL_TET_TRIPLETS){
            // checking one corner in this face suffices
            int i = ijk[0], j = ijk[1], k = ijk[2];
            if (smeared_curve_nc.scc.get(i, j, k) % 2 == 1){
                // there is an interior contractible region here; construct the interior region here
                // Lemma: can only exist as an interior hexagon
                if (smeared_curve_nc.scc.get(i, j, k) != 1 || smeared_curve_nc.scc.get(j, i, k) != 1 || smeared_curve_nc.scc.get(k, i, j) != 1){
                    throw std::logic_error("interior_region_construction: invalid contractible curve; should have exactly one segment on each face");
                }
                // ; the cv.k's here are helpers for mid scoop vertex later
                CombinatorialVertex cv_ij(i, j, 1, 2, k); // n_ij is 2 since there are two segments in this scenario
                CombinatorialVertex cv_ik(i, k, 1, 2, j);
                CombinatorialVertex cv_ji(j, i, 1, 2, k);
                CombinatorialVertex cv_jk(j, k, 1, 2, i);
                CombinatorialVertex cv_ki(k, i, 1, 2, j);
                CombinatorialVertex cv_kj(k, j, 1, 2, i);
                // Lemma: only one interior region can exist, so we can return after finding it
                return {CombFace({cv_ij, cv_ik, cv_ki, cv_kj, cv_jk, cv_ji}, CombFaceType::INTERIOR_HEXAGON)};
            }
        }
        return {}; // no interior region on any face
    }
    else {
        // Lemma: there can only be one interior region for a single corner-type curve
        int i = interior_vertex;
        // find the open side and emit the interior region from there
        // pick the open side (jk where scc.get(i,j,k) is 0) 
        int j, k; // choose the open jk
        int l;
        auto jkl = complement({i});
        if (smeared_curve_nc.scc.get(i, jkl[0], jkl[1]) == 0){
            j = jkl[0]; k = jkl[1];
        }
        else if (smeared_curve_nc.scc.get(i, jkl[0], jkl[2]) == 0){
            j = jkl[0]; k = jkl[2];
        }
        else if (smeared_curve_nc.scc.get(i, jkl[1], jkl[2]) == 0){
            j = jkl[1]; k = jkl[2];
        }
        else {
            throw std::invalid_argument("interior_region_construction: invalid corner curve; should have at least one open edge");
        }
        l = 6 - i - j - k; // the remaining vertex
        // construct the interior region; emit all scoop vertices on the jl edge; mid vertices will be added later
        // ; the cv.k's here are helpers for mid scoop vertex later
        // construct the relevant corner triangle as well; on corner i
        
        CombinatorialVertex cv_ij(i, j, 1, smeared_curve_nc.edge_ints[{i, j}], k);
        // go around the corner for the corner triangle
        CombinatorialVertex cv_il(i, l, 1, smeared_curve_nc.edge_ints[{i, l}], j);
        // get back on face ijk for the interior region
        CombinatorialVertex cv_ik(i, k, 1, smeared_curve_nc.edge_ints[{i, k}], j);
        CombinatorialVertex cv_kj(k, j, smeared_curve_nc.scc.get(k, j, i), smeared_curve_nc.edge_ints[{k, j}], i);
        vector<CombinatorialVertex> interior_region_vertices = {cv_ij, cv_il, cv_ik, cv_kj};
        // add scoop vertices between cv_kj and cv_jk
        int n_kj = smeared_curve_nc.edge_ints[{k, j}];
        for (int scoop_idx = 0; scoop_idx < smeared_curve_nc.scoop.get(i, j, k); scoop_idx++){
            CombinatorialVertex scoop_cv1(k, j, 2*scoop_idx + cv_kj.order + 1, n_kj, i); // cv_kj is first; and we got 1-based indexing
            CombinatorialVertex scoop_cv2(k, j, 2*scoop_idx + cv_kj.order + 2, n_kj, i);
            interior_region_vertices.push_back(scoop_cv1);
            interior_region_vertices.push_back(scoop_cv2);
        }
        CombinatorialVertex cv_jk(j, k, smeared_curve_nc.scc.get(j, k, i), smeared_curve_nc.edge_ints[{j, k}], i);
        interior_region_vertices.push_back(cv_jk);
        return {CombFace(interior_region_vertices, CombFaceType::INTERIOR_CORNER_TYPE)};
    }
}



// // CORNER TRIANGLE at the interior vertex (none for contractible curves).
// // Keeps scoop_mid_vertices because the vertex ordering depends on it.
// vector<CombFace>
// corner_triangle_construction(
//     const NormalCoordinates &smeared_curve_nc,
//     const int &interior_vertex,
//     bool scoop_mid_vertices
// ){
//     if (interior_vertex == -1){
//         return {}; // contractible curve; no corner triangle
//     }
//     int i = interior_vertex;
//     auto jkl = complement({i});
//     int j = jkl[0], k = jkl[1], l = jkl[2];

//     CombinatorialVertex cvij(i, j, 1, smeared_curve_nc.edge_ints[{i,j}]);
//     CombinatorialVertex cvik(i, k, 1, smeared_curve_nc.edge_ints[{i,k}]);
//     CombinatorialVertex cvil(i, l, 1, smeared_curve_nc.edge_ints[{i,l}]);
//     CombFace corner_comb_face;
//     if (!scoop_mid_vertices){
//         // just add the corner triangle; order doesnt matter
//         corner_comb_face = CombFace({cvij, cvik, cvil}, CombFaceType::CORNER_TYPE_TRIANGLE);
//     }
//     else {
//         // order matters here since mid vertex will be added on the interior end
//         // identify the open end where there is no corner segment; for potential scoop mid vertex addition
//         // for future construction: put the open part on the first edge of the triplet
//         if (smeared_curve_nc.scc.get(i, j, k) == 0){
//             corner_comb_face = CombFace({cvij, cvik, cvil}, CombFaceType::CORNER_TYPE_TRIANGLE);
//         }
//         else if (smeared_curve_nc.scc.get(i, j, l) == 0){
//             corner_comb_face = CombFace({cvij, cvil, cvik}, CombFaceType::CORNER_TYPE_TRIANGLE);
//         }
//         else if (smeared_curve_nc.scc.get(i, k, l) == 0){
//             corner_comb_face = CombFace({cvik, cvil, cvij}, CombFaceType::CORNER_TYPE_TRIANGLE);
//         }
//         else {
//             throw std::invalid_argument("corner_triangle_construction: invalid corner curve; should have at least one open edge");
//         }
//     }
//     // TODO: optionally add mid scoop vertex here as well? might look a bit nicer
//     return {corner_comb_face};
// }