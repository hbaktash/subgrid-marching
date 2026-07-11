#include "nc/NC_solver.h"


void 
NormalCoordinates::compute_NC(){

    // even-sum and triangle inequality checks done here
    compute_semi_corner_cuts_and_scoops();
    if ((even_sum && triangle_ineq) || !even_sum){ // if odd sum, just do the naive cut extraction
        SemiCornerCuts scc_reduced;
        compute_corner_cuts_and_reduced_semi_corner_cuts(scc_reduced);
        compute_diagonal_cuts(scc_reduced);
    }
    else if (even_sum){
        EdgeInts reduced_edge_ints = scoop_reduction();
        NormalCoordinates reduced_NC;
        reduced_NC.edge_ints = reduced_edge_ints;
        reduced_NC.compute_NC();
        if (!reduced_NC.triangle_ineq || !reduced_NC.even_sum){
            throw std::logic_error("Error: compute_NC: scoop reduction did not yield feasible NC!\n");
        }
        cc = reduced_NC.cc;
        dc = reduced_NC.dc;
    }
    else {
        return;
    }
}


void 
NormalCoordinates::compute_semi_corner_cuts_and_scoops(){
    for (std::array<int,3> triplet: ALL_TET_TRIPLETS){
        int i = triplet[0];
        int j = triplet[1];
        int k = triplet[2];
        int res_ij, res_ik, res_jk;
        bool feasible = triangle_inequality(edge_ints[{i,j}], edge_ints[{i,k}], edge_ints[{j,k}],
                                    res_ij, res_ik, res_jk);
        bool ijk_even_sum = ((edge_ints[{i,j}] + edge_ints[{i,k}] + edge_ints[{j,k}]) % 2 == 0);
        // log flags
        triangle_ineq = feasible && triangle_ineq;
        even_sum = ijk_even_sum && even_sum;
        
        if (!ijk_even_sum && feasible){
            res_ij = 1;
            res_ik = 1;
            res_jk = 1;
        }
        // Remove residuals/scoops
        int M_ij = edge_ints[{i,j}] - res_ij,
            M_ik = edge_ints[{i,k}] - res_ik,
            M_jk = edge_ints[{j,k}] - res_jk; 

        scc.set(i,j,k, std::max<int>( (M_ij + M_ik - M_jk) / 2, 0));
        scc.set(j,i,k, std::max<int>( (M_ij + M_jk - M_ik) / 2, 0));
        scc.set(k,i,j, std::max<int>( (M_ik + M_jk - M_ij) / 2, 0));

        // set scoops
        scoop.set(i,j,k, res_jk / 2);
        scoop.set(j,i,k, res_ik / 2);
        scoop.set(k,i,j, res_ij / 2);
        
    }
}

void
NormalCoordinates::compute_corner_cuts_and_reduced_semi_corner_cuts(SemiCornerCuts& scc_reduced){
    cc[0] = std::min({scc.get(0,1,2), scc.get(0,1,3), scc.get(0,2,3)});
    cc[1] = std::min({scc.get(1,0,2), scc.get(1,0,3), scc.get(1,2,3)});
    cc[2] = std::min({scc.get(2,0,1), scc.get(2,0,3), scc.get(2,1,3)});
    cc[3] = std::min({scc.get(3,0,1), scc.get(3,0,2), scc.get(3,1,2)});
    for (int i = 0; i < 4; i++){
        for (std::pair<int,int> jk_pair: all_pairs(complement({i}))){
            int j = jk_pair.first;
            int k = jk_pair.second;
            scc_reduced.set(i,j,k, scc.get(i,j,k) - cc[i]);
        }
    }
}

// reduced semicorner cuts to diagonal cuts
void
NormalCoordinates::compute_diagonal_cuts(const SemiCornerCuts& scc_reduced){
    // only iterate over {0,j} pairs
    for (int j = 1; j < 4; j++){
        std::vector<int> comp = complement({0,j});
        int k = comp[0];
        int l = comp[1];
        dc.set(0, j,
            std::min({
                    scc_reduced.get(k,0,j),
                    scc_reduced.get(l,0,j),
                    scc_reduced.get(j,k,l),
                    scc_reduced.get(0,k,l)
                }
            )
        );
    }
}

EdgeInts corner_and_diagonal_cuts_to_edge_ints(
    const CornerCuts& cc,
    const DiagonalCuts& dc
){
    EdgeInts edge_ints_ans;
    for (std::pair<int,int> ij_pair: all_pairs({0,1,2,3})){
        int i = ij_pair.first;
        int j = ij_pair.second;
        std::vector<int> kl_set = complement({i,j});
        int k = kl_set[0];
        int l = kl_set[1];
        edge_ints_ans[{i,j}] = cc[i] + cc[j] + dc[{i,k}] + dc[{i,l}];
    }
    return edge_ints_ans;
}

// get d1, d2
std::pair<int, int> NormalCoordinates::get_d1_d2() const {
    int b1 = 0;
    int diag_types = 0;
    int d1 = 0, d2 = 0;
    for (int i = 1; i < 4; i++){
        if (dc[{b1,i}] > 0) {
            diag_types++;
            if (d1 == 0) 
                d1 = dc[{b1,i}];
            else 
                d2 = dc[{b1,i}];
        }
    }
    if (diag_types > 2)
        throw std::logic_error("Error: get_d1_d2: more than 3 diagonal cut types, cannot determine d1, d2\n");
    return std::make_pair(d1, d2);
}

// more constructive stuff

EdgeInts 
NormalCoordinates::scoop_reduction(){
    // std::cout << "Scoop trace iteration...\n";
    EdgeOccupations edge_occupations(edge_ints); // sets all to false initially
    for (std::array<int,3> triplet: ALL_TET_TRIPLETS){
        int base, v1, v2;
        int i = triplet[0], j = triplet[1], k = triplet[2];
        if (scoop.get(i, j, k) > 0){
            base = i; v1 = j; v2 = k;
        }
        else if (scoop.get(j, k, i) > 0){
            base = j; v1 = k; v2 = i;
        }
        else if (scoop.get(k, i, j) > 0){
            base = k; v1 = i; v2 = j;
        }
        else {
            continue; // no scoops on this triangle
        }
        // generate combinatorial vertices for every scoop end follow it
        int scoop_count = scoop.get(base, v1, v2);
        for (int s = 0; s < scoop_count; s++){
            CombinatorialVertex cv_first(v1, v2, 2*s + scc.get(v1, v2, base) + 1, edge_ints[{v1, v2}], base);
            CombinatorialVertex cv_second(v1, v2, 2*s + scc.get(v1, v2, base) + 2, edge_ints[{v1, v2}], base);
            trace_forever(cv_first, edge_occupations);
            trace_forever(cv_second, edge_occupations);
        }
    }
    // make new edge ints from occupancies
    EdgeInts reduced_edge_ints;
    for (std::pair<int,int> ij_pair: all_pairs({0,1,2,3})){
        int i = ij_pair.first;
        int j = ij_pair.second;
        int occ_count = 0;
        for (size_t idx = 0; idx < edge_ints[{i,j}]; idx++){
            if (!edge_occupations.get(i,j,idx)){
                occ_count++;
            }
        }
        reduced_edge_ints.set(i,j, occ_count);
    }
    return reduced_edge_ints;
}

// helper to trace scoop curves
CombinatorialVertex 
NormalCoordinates::find_next_vertex(
    CombinatorialVertex v0
){
    if (even_sum == false)
        throw std::logic_error("combinatorial trace: even_sum is false, cannot trace scoops\n");
    int i = v0.i, j = v0.j;
    // new k on the other triangle
    int k = complement({i,j, v0.k})[0];
    // std::cout << "\t\t\t Tracing from vertex on edge {" << i << "," << j << "} at order " << v0.order << " with edge int " << v0.edge_int << "\n";
    // std::cout << "\t\t\t v0.k = " << v0.k << ", other vertex in face is " << k << "\n";

    int order = v0.order;
    int n_ij = v0.edge_int;

    // using i, j, k as the ordered vertices now
    if (scc.get(i, j, k) + scc.get(j, i, k)  + 2*scoop.get(k, i, j) != n_ij){
        throw std::logic_error("combinatorial trace: inconsistent scoop and semicorner cut data\n");
    }

    CombinatorialVertex v_next;
    if (order <= scc.get(i, j, k)){ // getting caught in i's corner turns
        // move to edge {i,k}
        v_next = CombinatorialVertex(i, k, order, edge_ints.get(i,k), j);
    }
    else if (order >= n_ij - scc.get(j, i, k) + 1){ // getting caught in j's corner turns
        // move to edge {j,k}
        v_next = CombinatorialVertex(j, k, n_ij - order + 1, edge_ints.get(j,k), i);
    }
    else{ // in scoop region
        // move to edge {i,j} but reverse direction
        int next_order;
        if ((order - scc.get(i, j, k)) % 2 == 0)
            next_order = order - 1;
        else 
            next_order = order + 1;
        v_next = CombinatorialVertex(i, j, next_order, n_ij, k);
    }
    return v_next;
}

void NormalCoordinates::trace_forever(CombinatorialVertex v0, EdgeOccupations &edge_occupations){
    CombinatorialVertex cv = v0;
    while (true){
        if (edge_occupations.get(cv)){
            return;
        }
        edge_occupations.set(cv, true);
        cv = find_next_vertex(cv);
    }
}


