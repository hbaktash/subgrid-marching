
#include "subgrid_MT/dual/dual_construction.h"
#include "subgrid_MT/dual/qef_solver.h"
#include "subgrid_MT/cv_interpolation.h"

using std::vector;


TriangleSoup
dual_subgrid_surface(
    const std::array<Vector3,4>& tet_positions,
    const std::array<size_t,4>& tet_indices,
    const std::array<std::vector<double>,6>& edge_isect_ts,
    const std::array<std::vector<Vector3>,6>& edge_isect_normals,
    const std::vector<CombFace>& open_curves,
    const std::vector<CombFace>& scoop_curves,
    const std::vector<CombFace>& normal_curves,
    double reg_alpha,
    bool project_duals
){
    vector<CombFace> curves;
    for (const auto& c : scoop_curves) curves.push_back(c);
    for (const auto& c : normal_curves) curves.push_back(c);
    if (curves.empty()) return {};

    auto normals = normals_from_cv_polygons(curves, edge_isect_ts, edge_isect_normals);
    TriangleSoup soup = build_dual_local_soup(curves, tet_positions, edge_isect_ts);

    soup.dual_positions = dual_positions_from_isect_data(
        tet_positions, soup.faces, soup.vertices, normals,
        reg_alpha, project_duals
    );

    std::array<int,4> idx_map = {(int)tet_indices[0], (int)tet_indices[1],
                                  (int)tet_indices[2], (int)tet_indices[3]};
    soup.remap_vertex_indices(idx_map);
    return soup;
}


std::vector<Vector3>
normals_from_cv_polygons(
    const std::vector<CombFace>& curves,
    const std::array<std::vector<double>,6>& edge_isect_ts,
    const std::array<std::vector<Vector3>,6>& edge_isect_normals
){
    std::vector<Vector3> normals;
    for (const auto& curve : curves){
        for (const auto& cv : curve.vertices){
            int i = cv.i, j = cv.j;
            int edge_idx = edge_pair_to_index(i, j);
            size_t idx = (i < j) ? (cv.order - 1) : (cv.edge_int - cv.order);
            normals.push_back(edge_isect_normals[edge_idx][idx]);
        }
    }
    return normals;
}


TriangleSoup
build_dual_local_soup(
    const std::vector<CombFace>& curves,
    const std::array<Vector3,4>& tet_positions,
    const std::array<std::vector<double>,6>& edge_isect_ts
){
    TriangleSoup soup;
    size_t vert_offset = 0;
    for (const auto& curve : curves){
        size_t n = curve.vertices.size();
        vector<size_t> polygon(n);
        auto fpe = curve_segment_faces(curve);
        for (size_t vi = 0; vi < n; vi++){
            const auto& cv = curve.vertices[vi];
            soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0));
            soup.signatures.push_back(key_from_cv(cv));
            polygon[vi] = vert_offset + vi;
        }
        soup.faces.push_back(polygon);
        soup.faces_per_edge.push_back(fpe);
        vert_offset += n;
    }
    return soup;
}
