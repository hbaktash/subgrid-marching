#include "subgrid_MT/normal_reconstruction.h"

#include <numeric>


// Construction policy (NOT a runtime option): separate the corner disks from the
// recursive subdivision. This is the shipped behavior and should only be changed
// by a developer who understands the almost-normal recursion. It is deliberately
// a file-local constant rather than a plumbed-through flag.
static constexpr bool no_corner_recursion = true;


static void shift_orders_past_corners(
    TriangleSoup &soup,
    const NormalCoordinates &NC
){
    for (auto &sig : soup.signatures){
        if (sig.type == CombVertexSigType::STANDALONE) continue;
        if (sig.i == 4 || sig.j == 4) continue;
        sig.order += NC.cc[sig.i];
        sig.edge_int = NC.edge_ints[{sig.i, sig.j}];
    }
}


TriangleSoup
normal_surface_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    bool verbose
){
    EdgeInts normal_edge_ints(normal_edge_isect_ts);
    NormalCoordinates NC(normal_edge_ints);
    if (!NC.triangle_ineq || !NC.even_sum){
        // print the edge ints
        if (verbose){
            std::cout << "Invalid normal coordinates: triangle inequality or even sum condition failed." << std::endl;
            normal_edge_ints.print();
        }
        throw std::runtime_error("Invalid normal coordinates: triangle inequality or even sum condition failed.");
    }

    auto [d1, d2] = NC.get_d1_d2();

    if (d1 != 0 && d2 != 0){      if (verbose)      std::cout << "Almost-normal Surface ---- d1 = " << d1 << ", d2 = " << d2 << std::endl;}
    else if (d1 != 0 || d2 != 0){ if (verbose) std::cout << "Normal Surface ---- d1 = " << d1 << ", d2 = " << d2 << std::endl;}
    else if (d1 == 0 && d2 == 0){ if (verbose) std::cout << "No diags ---- d1 == d2 == 0\n"; }
    else{ throw std::runtime_error("normal_boundary_surface: Unexpected case for d1, d2 values; make sure intersections are normal."); }

    TriangleSoup soup;
    if (d1 == 0 || d2 == 0){ // Normal surface; triangles, quads
        auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
        auto diagonal_soup = diagonal_construction(tet_positions, normal_edge_isect_ts, NC);
        soup.add_local_soup(corner_soup);
        soup.add_local_soup(diagonal_soup);
    }
    else { // Almost-normal surface; octagons, spirals, recursion
        if (d1 == d2 || std::gcd(d1, d2) == 1){ // no recursion
            auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
            soup.add_local_soup(corner_soup);
            if (d1 == d2){
                auto octagon_soup = octagon_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(octagon_soup);
            }
            else {
                auto spiral_soup = spiral_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(spiral_soup);
            }
        }
        else { // recursion
            if (verbose) std::cout << "Recursion ---- d1 = " << d1 << ", d2 = " << d2 << std::endl;
            if (no_corner_recursion){
                auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(corner_soup);
            }
            auto recursive_soup = recursive_construction(tet_positions, normal_edge_isect_ts, NC, verbose);
            soup.add_local_soup(recursive_soup);
        }
    }
    return soup;
}



// ------------------ Normal Surface stuff ------------------
TriangleSoup
corner_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
){
    TriangleSoup soup;
    for (int i = 0; i < 4; i++){
        std::vector<int> others = complement({i});
        int j = others[0];
        int k = others[1];
        int l = others[2];
        for (int order = 1; order <= NC.cc[i]; order++){
            vector<CombinatorialVertex> triangle_comb_verts{
                CombinatorialVertex(i, j, order, NC.edge_ints[{i,j}]),
                CombinatorialVertex(i, k, order, NC.edge_ints[{i,k}]),
                CombinatorialVertex(i, l, order, NC.edge_ints[{i,l}])
            };
            CombFace comb_face(triangle_comb_verts, CombFaceType::TRIANGLE);
            soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
        }
    }
    return soup;
}

TriangleSoup
diagonal_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
){
    TriangleSoup soup;
    int b1 = 0;
    int b2 = NC.dc[{b1, 1}] != 0 ? 1 : NC.dc[{b1, 2}] != 0 ? 2 : NC.dc[{b1, 3}] != 0 ? 3 : -1;
    if (b2 == -1) return soup; // no diagonal cuts

    std::vector<int> head_vs = complement({b1, b2});
    int h1 = head_vs[0], h2 = head_vs[1];

    int diag_count = NC.dc[{b1, b2}];
    int b1_cc = NC.cc[b1],
        b2_cc = NC.cc[b2];

    for(int t = 0; t < diag_count; t++){
        CombinatorialVertex cv0(b1, h1, t + b1_cc + 1, NC.edge_ints[{b1, h1}]);
        CombinatorialVertex cv1(b1, h2, t + b1_cc + 1, NC.edge_ints[{b1, h2}]);
        CombinatorialVertex cv2(b2, h2, t + b2_cc + 1, NC.edge_ints[{b2, h2}]);
        CombinatorialVertex cv3(b2, h1, t + b2_cc + 1, NC.edge_ints[{b2, h1}]);
        CombFace comb_face({cv0, cv1, cv2, cv3}, CombFaceType::QUAD02);
        soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    }
    return soup;
}



// ------------------ Almost Normal Surface stuff ------------------
TriangleSoup
octagon_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
){
    TriangleSoup soup;
    int b1 = 0;
    std::vector<int> other_vs{1, 2, 3};
    int vd1 = -1, vd2 = -1, vdd = -1;
    int d1 = -1, d2 = -1;
    for(int other_b1: other_vs){
        int tmp_d = NC.dc[{b1, other_b1}];
        if(tmp_d != 0){
            if(vd1 == -1){
                vd1 = other_b1;
                d1 = tmp_d;
                continue;
            }
            else {
                vd2 = other_b1;
                d2 = tmp_d;
            }
        }
        else vdd = other_b1;
    }
    assert(vd1 != -1 && vd2 != -1 && vdd != -1 && d1 > 0 && d2 > 0 && d1 == d2);

    int d = d1;
    int c_b1 = NC.cc[b1], c_vd1 = NC.cc[vd1], c_vd2 = NC.cc[vd2], c_vdd = NC.cc[vdd];

    CombinatorialVertex b1_vdd_mid_vert_1(b1, vdd, c_b1 + d, NC.edge_ints[{b1, vdd}]);
    CombinatorialVertex b1_vdd_mid_vert_2(b1, vdd, c_b1 + d + 1, NC.edge_ints[{b1, vdd}]);
    DerivedVertex b1_vdd_mid_vertex({b1_vdd_mid_vert_1, b1_vdd_mid_vert_2});

    CombinatorialVertex vd1_vd2_mid_vert_1(vd1, vd2, c_vd1 + d, NC.edge_ints[{vd1, vd2}]);
    CombinatorialVertex vd1_vd2_mid_vert_2(vd1, vd2, c_vd1 + d + 1, NC.edge_ints[{vd1, vd2}]);
    DerivedVertex vd1_vd2_mid_vertex({vd1_vd2_mid_vert_1, vd1_vd2_mid_vert_2});

    for(int t = 0; t < d; t++){
        DerivedVertex inner_vertex_t(vd1_vd2_mid_vertex, b1_vdd_mid_vertex, (double)(t + 1)/(double)(d + 1));

        CombinatorialVertex b1vdd_1(b1 , vdd, t + 1 + c_b1  , NC.edge_ints[{b1, vdd}] );
        CombinatorialVertex b1vdd_2(vdd, b1 , t + 1 + c_vdd , NC.edge_ints[{b1, vdd}] );
        CombinatorialVertex d1d2_1 (vd1, vd2, d - t + c_vd1 , NC.edge_ints[{vd1, vd2}]);
        CombinatorialVertex d1d2_2 (vd2, vd1, d - t + c_vd2 , NC.edge_ints[{vd1, vd2}]);
        CombinatorialVertex b1vd1  (b1 , vd1, t + 1 + c_b1  , NC.edge_ints[{b1, vd1}] );
        CombinatorialVertex b1vd2  (b1 , vd2, t + 1 + c_b1  , NC.edge_ints[{b1, vd2}] );
        CombinatorialVertex d1vdd  (vd1, vdd, d - t + c_vd1 , NC.edge_ints[{vd1, vdd}]);
        CombinatorialVertex d2vdd  (vd2, vdd, d - t + c_vd2 , NC.edge_ints[{vd2, vdd}]);

        CombFace octagon_comb_face({
            b1vdd_1, b1vd1, d1d2_1, d1vdd, b1vdd_2, d2vdd, d1d2_2, b1vd2
        }, CombFaceType::OCTAGON);
        octagon_comb_face.derived_vertex = inner_vertex_t;

        soup.add_local_soup(embed_comb_face(octagon_comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    }
    return soup;
}


TriangleSoup spiral_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
){
    assert(NC.even_sum && NC.triangle_ineq);

    int i = 0, j = 1;
    CombinatorialVertex start_cv(i, j, NC.cc[i] + 1, NC.edge_ints[{i,j}]);

    vector<CombinatorialVertex> spiral_comb_verts;
    std::array<size_t, 4> dummy_global_inds = {0, 1, 2, 3};
    EdgeOccupations dummy_edge_occ;
    trace_generic_curve(start_cv, NC, dummy_global_inds, dummy_edge_occ, spiral_comb_verts, false);

    CombFace spiral_comb_face(spiral_comb_verts, CombFaceType::SPIRAL);
    TriangleSoup soup;
    soup.add_local_soup(embed_comb_face(spiral_comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    return soup;
}



TriangleSoup recursive_construction(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC,
    bool verbose
){
    auto [d1, d2] = NC.get_d1_d2();
    if (d1 < d2) std::swap(d1, d2);
    std::array<int, 4> interior_edge_ints = {0, 0, 0, 0};
    int b1 = 0;
    interior_edge_ints[b1] = d1 - d2 + (no_corner_recursion ? 0 : NC.cc[b1]);
    std::array<int, 3> others_vs = {1, 2, 3};
    for (int other_v : others_vs){
        if (NC.dc[{b1, other_v}] == 0) interior_edge_ints[other_v] = 2*d2 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else if (NC.dc[{b1, other_v}] == d1) interior_edge_ints[other_v] = d1 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else if (NC.dc[{b1, other_v}] == d2) interior_edge_ints[other_v] = d2 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else throw std::logic_error("recursive_construction: diagonal cuts cant take any other values. NC values are off\n");
    }

    Vector3 new_interior_tet_vertex = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0;
    std::array<double, 4> interior_subd_ts_limit = {1.0, 1.0, 1.0, 1.0};
    
    if (no_corner_recursion){
        auto [recursion_center, interior_isect_limits] = compute_recursion_center_and_isect_ts(tet_positions, normal_edge_isect_ts, NC);
        new_interior_tet_vertex = recursion_center;
        for (int i = 0; i < 4; i++){
            interior_subd_ts_limit[i] = interior_isect_limits[i];
        }
    }
        

    TriangleSoup soup;
    for (std::array<int, 3> triplet: ALL_TET_TRIPLETS){
        int i = triplet[0], j = triplet[1], k = triplet[2];

        array<Vector3,4> interior_tet_positions = {
            new_interior_tet_vertex, tet_positions[i], tet_positions[j], tet_positions[k]
        };

        array<vector<double>, 6> interior_tet_edge_isect_ts;
        interior_tet_edge_isect_ts[0] = generate_uniform_isect_ts_single_edge(interior_edge_ints[i], interior_subd_ts_limit[i]);
        interior_tet_edge_isect_ts[1] = generate_uniform_isect_ts_single_edge(interior_edge_ints[j], interior_subd_ts_limit[j]);
        interior_tet_edge_isect_ts[2] = generate_uniform_isect_ts_single_edge(interior_edge_ints[k], interior_subd_ts_limit[k]);
        if (no_corner_recursion){
            auto& ts_ij = normal_edge_isect_ts[edge_pair_to_index(i, j)];
            auto& ts_ik = normal_edge_isect_ts[edge_pair_to_index(i, k)];
            auto& ts_jk = normal_edge_isect_ts[edge_pair_to_index(j, k)];
            interior_tet_edge_isect_ts[3] = vector<double>(ts_ij.begin() + NC.cc[i], ts_ij.end() - NC.cc[j]);
            interior_tet_edge_isect_ts[4] = vector<double>(ts_ik.begin() + NC.cc[i], ts_ik.end() - NC.cc[k]);
            interior_tet_edge_isect_ts[5] = vector<double>(ts_jk.begin() + NC.cc[j], ts_jk.end() - NC.cc[k]);
        } else {
            interior_tet_edge_isect_ts[3] = normal_edge_isect_ts[edge_pair_to_index(i, j)]; // this order is fine since ALL_TET_TRIPLETS returns i < j < k
            interior_tet_edge_isect_ts[4] = normal_edge_isect_ts[edge_pair_to_index(i, k)];
            interior_tet_edge_isect_ts[5] = normal_edge_isect_ts[edge_pair_to_index(j, k)];
        }

        auto sub_soup = normal_surface_construction(interior_tet_positions, interior_tet_edge_isect_ts, verbose);
        if (TriangleSoup::COMB_MERGE){
            sub_soup.remap_vertex_indices({4, i, j, k});
            if (no_corner_recursion)
                shift_orders_past_corners(sub_soup, NC);
        }
        soup.add_local_soup(sub_soup);
    }
    if (TriangleSoup::COMB_MERGE)
        soup.make_vertex_standalone(4);
    return soup;
}




// function for determining the recursion center and interior isect_ts, when no_corner_recursion is true
std::pair<Vector3, std::array<double, 4>> 
compute_recursion_center_and_isect_ts(
    const std::array<Vector3,4> &tet_positions,
    const std::array<std::vector<double>,6> &normal_edge_isect_ts,
    const NormalCoordinates &NC
){
    // get the first non-corner cv's at each edge; from either side
    CombinatorialVertex cv01(0, 1, NC.cc[0] + 1, NC.edge_ints[{0, 1}]);
    CombinatorialVertex cv10(1, 0, NC.cc[1] + 1, NC.edge_ints[{0, 1}]);
    CombinatorialVertex cv02(0, 2, NC.cc[0] + 1, NC.edge_ints[{0, 2}]);
    CombinatorialVertex cv20(2, 0, NC.cc[2] + 1, NC.edge_ints[{0, 2}]);
    CombinatorialVertex cv03(0, 3, NC.cc[0] + 1, NC.edge_ints[{0, 3}]);
    CombinatorialVertex cv30(3, 0, NC.cc[3] + 1, NC.edge_ints[{0, 3}]);
    CombinatorialVertex cv12(1, 2, NC.cc[1] + 1, NC.edge_ints[{1, 2}]);
    CombinatorialVertex cv21(2, 1, NC.cc[2] + 1, NC.edge_ints[{1, 2}]);
    CombinatorialVertex cv13(1, 3, NC.cc[1] + 1, NC.edge_ints[{1, 3}]);
    CombinatorialVertex cv31(3, 1, NC.cc[3] + 1, NC.edge_ints[{1, 3}]);
    CombinatorialVertex cv23(2, 3, NC.cc[2] + 1, NC.edge_ints[{2, 3}]);
    CombinatorialVertex cv32(3, 2, NC.cc[3] + 1, NC.edge_ints[{2, 3}]);
    // convert to positions
    Vector3 p01 = interpolate_comb_vertex(cv01, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p10 = interpolate_comb_vertex(cv10, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p02 = interpolate_comb_vertex(cv02, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p20 = interpolate_comb_vertex(cv20, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p03 = interpolate_comb_vertex(cv03, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p30 = interpolate_comb_vertex(cv30, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p12 = interpolate_comb_vertex(cv12, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p21 = interpolate_comb_vertex(cv21, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p13 = interpolate_comb_vertex(cv13, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p31 = interpolate_comb_vertex(cv31, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p23 = interpolate_comb_vertex(cv23, tet_positions, normal_edge_isect_ts, -1.0);
    Vector3 p32 = interpolate_comb_vertex(cv32, tet_positions, normal_edge_isect_ts, -1.0);

    // subdivision center is the average of these 12 points
    Vector3 subd_center = (p01 + p10 + p02 + p20 + p03 + p30 + p12 + p21 + p13 + p31 + p23 + p32) / 12.0;

    // now determine the t, where the segment from the subd_center to each tet vertex intersects the corresponding corner triangle formed by thse 12 points
    // segment 0: subd_center -> tet_positions[0]; triangle: p01, p02, p03
    double ts0_limit, ts1_limit, ts2_limit, ts3_limit;
    bool hit0 = segment_triangle_intersection(subd_center, tet_positions[0], p01, p02, p03, ts0_limit);
    bool hit1 = segment_triangle_intersection(subd_center, tet_positions[1], p10, p12, p13, ts1_limit);
    bool hit2 = segment_triangle_intersection(subd_center, tet_positions[2], p20, p21, p23, ts2_limit);
    bool hit3 = segment_triangle_intersection(subd_center, tet_positions[3], p30, p31, p32, ts3_limit);
    if (!hit0 || !hit1 || !hit2 || !hit3){
        log_warn("compute_recursion_center_and_isect_ts: segment-triangle intersection failed; falling back to greedy recursion intersections.");
        // Fallback to greedy recursion intersections
        ts0_limit = ts1_limit = ts2_limit = ts3_limit = 1.0;
    }

    return std::make_pair(subd_center, std::array<double, 4>{{ts0_limit, ts1_limit, ts2_limit, ts3_limit}});
}