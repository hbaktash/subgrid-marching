#include "subgrid_MT/primal_reconstruction.h"
#include "subgrid_MT/simplicial_embedding.h"

// picks a subset of intersections to only unoccupied ones
std::array<std::vector<double>,6> 
unoccupied_edge_isect_ts(
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const EdgeOccupations &edge_occupations
){
    std::array<std::vector<double>,6> result;
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        int i = ALL_TET_PAIRS[eIdx].first,
            j = ALL_TET_PAIRS[eIdx].second;
        int n_ij = edge_isect_ts[eIdx].size();
        for (int k = 0; k < n_ij; ++k) {
            if (!edge_occupations.get(i, j, k)) {
                result[eIdx].push_back(edge_isect_ts[eIdx][k]);
            }
        }
    }
    return result;
}

TriangleSoup
subgrid_greedy(
    const std::array<Vector3,4> &tet_positions,
    const std::array<size_t,4> &tet_global_indices,
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const vector<CombFace> &scoop_boundary_curves,
    const vector<CombFace> &normal_boundary_curves,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge
){
    TriangleSoup soup;

    auto embed_greedy = [&](CombFace curve) {
        if (curve.vertices.size() == 3) {
            curve.type = CombFaceType::TRIANGLE;
        } else {
            curve.type = CombFaceType::SPIRAL;
        }
        soup.add_local_soup(embed_comb_face(curve, tet_positions, edge_isect_ts, false, scoop_mid_vertex_bulge));
    };

    for (const auto &curve : normal_boundary_curves)
        embed_greedy(curve);

    for (auto curve : scoop_boundary_curves) {
        if (scoop_mid_vertices)
            curve = add_mid_scoop_vertices_spiral(curve);
        embed_greedy(curve);
    }

    if (TriangleSoup::COMB_MERGE) {
        array<int,4> idx_map = {(int)tet_global_indices[0], (int)tet_global_indices[1],
                                 (int)tet_global_indices[2], (int)tet_global_indices[3]};
        soup.remap_vertex_indices(idx_map);
        soup.make_type_standalone(CombVertexSigType::SCOOP_INTERIOR_STEINER);
    }
    return soup;
}


TriangleSoup
subgrid_surface(
    const std::array<Vector3,4> &tet_positions,
    const std::array<size_t,4> &tet_global_indices,
    const std::array<std::vector<double>,6> &edge_isect_ts,
    const vector<CombFace> &open_boundary_curves,
    const vector<CombFace> &scoop_boundary_curves,
    const EdgeOccupations &non_normal_edge_occupations,
    bool scoop_mid_vertices,
    double scoop_mid_vertex_bulge,
    bool verbose
){
    TriangleSoup soup;
    // open curves are discarded
    
    // evensum surface construction
    TriangleSoup evensum_soup = evensum_surface_construction(tet_positions, edge_isect_ts, scoop_boundary_curves,
                scoop_mid_vertices, scoop_mid_vertex_bulge);
    soup.add_local_soup(evensum_soup);

    // normal surface construction
    auto normal_edge_isect_ts = unoccupied_edge_isect_ts(edge_isect_ts, non_normal_edge_occupations);
    auto normal_soup = normal_surface_construction(tet_positions, normal_edge_isect_ts, verbose);
    if (TriangleSoup::COMB_MERGE)
        normal_soup.remap_orders(edge_isect_ts, non_normal_edge_occupations, false);
    soup.add_local_soup(normal_soup);
    if (TriangleSoup::COMB_MERGE){
        array<int,4> idx_map = {(int)tet_global_indices[0], (int)tet_global_indices[1],
                                 (int)tet_global_indices[2], (int)tet_global_indices[3]};
        soup.remap_vertex_indices(idx_map);
        soup.make_type_standalone(CombVertexSigType::SCOOP_INTERIOR_STEINER);
    }
    return soup;
}

