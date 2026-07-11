#include "subgrid_MT/boundary_curve.h"

bool
pick_residual_cv(
    int i, int j, int other,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds,
    CombinatorialVertex &out_cv
){
    int nij = nc.edge_ints[{i,j}];
    int scoop_ij = nc.scoop.get(other, i, j);
    int scc_i = nc.scc.get(i,j,other);
    int scc_j = nc.scc.get(j,i,other);

    if (nij < scoop_ij * 2 + scc_i + scc_j){
        throw std::logic_error("pick_residual_cv: inconsistent scoop and scc counts\n");
    }
    else if (nij == scoop_ij * 2 + scc_i + scc_j){
        return false; // no residual cv on this edge
    }
    else if (nij == scoop_ij * 2 + scc_i + scc_j + 1){ // expect exactly 1 residual cv on any edge
        if (scoop_ij == 0){
            // just pick the non-corner vertex
            int order = nc.scc.get(i,j,other) + 1;
            out_cv = CombinatorialVertex(i, j, order, nij, other);
        }
        else { // must pick in the scoop region
            // even is easy, just pick the middle
            if (scoop_ij % 2 == 0){
                int order = nc.scc.get(i,j,other) + scoop_ij + 1; // scoop_ij is half of scoop vertices
                out_cv = CombinatorialVertex(i, j, order, nij, other);
            }
            else { // odd scoop count, pick the vertex closer to smaller index
                int order = nc.scc.get(i,j,other) + scoop_ij - 1 + 1;
                if (global_inds[i] > global_inds[j])
                    order += 2;
                out_cv = CombinatorialVertex(i, j, order, nij, other);
            }
        }
        return true;
    }
    else {
        throw std::logic_error("pick_residual_cv: more than 1 residual cv on this edge, cannot pick uniquely\n");
    }
}


tuple<vector<CombFace>, vector<CombFace>, vector<CombFace>>
boundary_comb_curves(
    const array<size_t,4> tet_global_indices,
    const array<int, 6> edge_isect_counts,
    EdgeOccupations &out_edge_occupations,
    bool populate_normal_curves
){
    // pre compute 2D/semi-corner and scoop data
    EdgeInts edge_ints{edge_isect_counts};
    NormalCoordinates NC;
    NC.edge_ints = edge_ints;
    NC.compute_semi_corner_cuts_and_scoops(); // only need scoops and semi corner data
    
    // pick starting CV's; residuals (odd), scoops (even)
    vector<CombinatorialVertex> residual_cvs, scoop_cvs;
    pick_residuals_and_scoop_cvs(NC, tet_global_indices, residual_cvs, scoop_cvs);

    // Trace all curves 
    // start from odd curves; since they are open it matters to start from an open end.
    EdgeOccupations edge_occupations(edge_ints); // defaults to all false
    vector<CombFace> open_curves;
    for (const auto &res_cv: residual_cvs){
        if (!edge_occupations.get(res_cv)){ // could be the tail of an already traced curve
            vector<CombinatorialVertex> curve_vertices;
            trace_generic_curve(res_cv, NC, tet_global_indices, edge_occupations, curve_vertices);
            open_curves.emplace_back(curve_vertices, CombFaceType::OPEN);
        }
    }
    
    // trace scooped/even-sum curves
    vector<CombFace> scoop_curves;
    for (const auto &scoop_cv: scoop_cvs){
        if (!edge_occupations.get(scoop_cv)){ // could be part of an already traced OPEN curve; TODO: actually might not be possible, so check for future work reasons
            vector<CombinatorialVertex> curve_vertices;
            trace_generic_curve(scoop_cv, NC, tet_global_indices, edge_occupations, curve_vertices);
            scoop_curves.emplace_back(curve_vertices, CombFaceType::CLOSED_NON_NORMAL);
        }
    }

    // return for primal normal construction
    out_edge_occupations = edge_occupations; 

    // normal curves; not needed for primal construction; only for dual and visuals..
    vector<CombFace> normal_curves;
    if (populate_normal_curves){
        for (auto ij: ALL_TET_PAIRS){
            int i = ij.first, j = ij.second;
            for (int order = 1; order <= edge_ints[{i,j}]; order++){ // inefficient, but ok
                CombinatorialVertex cv(i, j, order, edge_ints[{i,j}]); // convert to 1-based order
                if (!edge_occupations.get(cv)){
                    vector<CombinatorialVertex> curve_vertices;
                    trace_generic_curve(cv, NC, tet_global_indices, edge_occupations, curve_vertices);
                    normal_curves.emplace_back(curve_vertices, CombFaceType::CLOSED_NORMAL);
                }
            }
        }
    }
    return {open_curves, scoop_curves, normal_curves};
}


void 
trace_generic_curve(
    const CombinatorialVertex &start_cv,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds,
    EdgeOccupations &occupations,
    std::vector<CombinatorialVertex> &curve_vertices,
    bool keep_occupation_track // = true
){
    // Caller should check occupation of start_cv before calling this function
    CombinatorialVertex current_cv = start_cv;
    if (keep_occupation_track) occupations.set(current_cv, true);
    curve_vertices.push_back(current_cv);
    // std::cout << " \t inside trace_generic_curve\n";
    while (true){
        CombinatorialVertex next_cv = find_next_generic_vertex(current_cv, nc, global_inds);
        // std::cout << "next vertex: \n\t";
        // next_cv.print();
        // std::cout << "\n";
        if (keep_occupation_track) occupations.set(next_cv, true);
        // check termination
        if (next_cv == current_cv || next_cv == start_cv){
            break;
        }
        // add to curve
        curve_vertices.push_back(next_cv);
        current_cv = next_cv;
    }
}



CombinatorialVertex 
find_next_generic_vertex(
    const CombinatorialVertex &cv0,
    const NormalCoordinates &nc,
    const std::array<size_t, 4> &global_inds
){
    int i = cv0.i, j = cv0.j;
    // new k on the other triangle
    int k = complement({i,j, cv0.k})[0];
    int order = cv0.order;
    int n_ij = cv0.edge_int;
    CombinatorialVertex v_next;

    // if its on the corner cuts, then follow like usual
    if (order <= nc.scc.get(i, j, k)){ // getting caught in i's corner turns
        // move to edge {i,k}
        v_next = CombinatorialVertex(i, k, order, nc.edge_ints.get(i,k), j);
        return v_next;
    }
    else if (order >= n_ij - nc.scc.get(j, i, k) + 1){ // getting caught in j's corner turns
        // move to edge {j,k}
        v_next = CombinatorialVertex(j, k, n_ij - order + 1, nc.edge_ints.get(j,k), i);
        return v_next;
    }

    // otherwise, in scoop region, need to check for residual vertex
    CombinatorialVertex cv_residual;
    bool has_residual = pick_residual_cv(i, j, k, nc, global_inds, cv_residual);
    int res_order = has_residual ? cv_residual.order : -1;

    if (order == res_order){
        // terminated at residual vertex
        return cv_residual;
    }
    else {
        int next_order;
        if (order < res_order || !has_residual){ // either no residual or passed it; so no effect
            if ((order - nc.scc.get(i, j, k)) % 2 == 0)
                next_order = order - 1;
            else 
                next_order = order + 1;
        }
        else {
            if ((order - nc.scc.get(i, j, k)) % 2 == 0)
                next_order = order + 1;
            else 
                next_order = order - 1;
        }
        v_next = CombinatorialVertex(i, j, next_order, n_ij, k);
        return v_next;
    }
}



void pick_residuals_and_scoop_cvs(
    const NormalCoordinates &NC,
    const array<size_t,4> &tet_global_indices,
    vector<CombinatorialVertex> &residual_cvs,
    vector<CombinatorialVertex> &scoop_cvs
){
    for (const auto &ijk: ALL_TET_TRIPLETS) {
        for (auto ij: all_pairs(ijk)){
            int i = ij.first, j = ij.second;
            int k = ijk[0] + ijk[1] + ijk[2] - i - j;
            CombinatorialVertex residual_cv;
            if (pick_residual_cv(i, j, k, NC, tet_global_indices, residual_cv)){
                residual_cvs.push_back(residual_cv);
            }
            auto ij_scoop_cvs = pick_scoop_cvs(i, j, k, NC, residual_cv);
            scoop_cvs.insert(scoop_cvs.end(), ij_scoop_cvs.begin(), ij_scoop_cvs.end());
        }
    }
}


std::vector<CombinatorialVertex>
pick_scoop_cvs(
    int i, int j, int other,
    const NormalCoordinates &nc,
    const CombinatorialVertex &residual_cv
){
    assert((nc.edge_ints[{i,j}] >= 2 * nc.scoop.get(other, i, j) + nc.scc.get(i,j,other) + nc.scc.get(j,i,other)));
    
    int nij = nc.edge_ints[{i,j}];
    bool has_residual = residual_cv.i == -1; // default cv 
    
    // go over all possible vertices on this edge from one scc to the other
    int start_order = nc.scc.get(i,j,other) + 1,
        end_order = nij - nc.scc.get(j,i,other);
    vector<CombinatorialVertex> scoop_cvs;

    for (int order = start_order; order <= end_order; order++){
        CombinatorialVertex cv(i, j, order, nij, other);
        if (!has_residual || !(cv == residual_cv)){
            scoop_cvs.push_back(cv);
        }
    }
    return scoop_cvs;
}


std::array<size_t,3>
segment_tet_face(
    const CombinatorialVertex &v_a,
    const CombinatorialVertex &v_b
){
    if (v_a.same_edge_as(v_b)) {
        return {(size_t)v_b.i, (size_t)v_b.j, (size_t)v_b.k};
    }
    int third = (v_a.i != v_b.i && v_a.i != v_b.j) ? v_a.i : v_a.j;
    return {(size_t)v_b.i, (size_t)v_b.j, (size_t)third};
}


std::vector<std::array<size_t,3>>
curve_segment_faces(
    const CombFace &comb_face
){
    const auto &verts = comb_face.vertices;
    bool is_closed = (comb_face.type != CombFaceType::OPEN);
    size_t n = verts.size();
    size_t n_segments = is_closed ? n : n - 1;
    std::vector<std::array<size_t,3>> faces;
    faces.reserve(n_segments);
    for (size_t i = 0; i < n_segments; ++i)
        faces.push_back(segment_tet_face(verts[i], verts[(i + 1) % n]));
    return faces;
}
