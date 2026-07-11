#include "subgrid_MT/simplicial_embedding.h"


void
align_cvs(
    const CombinatorialVertex &cv1,
    const CombinatorialVertex &cv2,
    CombinatorialVertex &aligned_cv1,
    CombinatorialVertex &aligned_cv2
){
    if (!cv1.same_edge_as(cv2))
        throw std::invalid_argument("align_cvs: vertices are not on the same edge");
    if (cv1.edge_int != cv2.edge_int)
        throw std::invalid_argument("align_cvs: vertices do not have the same edge intersection count");
    // align the combinatorial vertices
    if (cv1.i == cv2.i && cv1.j == cv2.j){
        // same direction; align orders
        aligned_cv1 = cv1;
        aligned_cv2 = cv2;
    }
    else if (cv1.i == cv2.j && cv1.j == cv2.i){
        int n_ij = cv1.edge_int; // same as cv2.edge_int
        // opposite direction; align orders and swap i,j
        aligned_cv1 = CombinatorialVertex(cv1.i, cv1.j, cv1.order, n_ij, cv1.k);
        aligned_cv2 = CombinatorialVertex(cv1.i, cv1.j, n_ij - cv2.order + 1, n_ij, cv2.k);
    }
    else {
        throw std::logic_error("align_cvs: invalid vertices");
    }
}


vector<int>
get_close_scoop_sides(
    const NormalCoordinates &smeared_curve_nc,
    const CombinatorialVertex &inner_segment_cv1,
    const CombinatorialVertex &inner_segment_cv2
){
    if(!inner_segment_cv1.same_edge_as(inner_segment_cv2)){
        throw std::invalid_argument("get_close_scoop_sides: vertices are not on the same edge");
    }
    int i = inner_segment_cv1.i, j = inner_segment_cv1.j,
        n_ij = smeared_curve_nc.edge_ints[{i, j}];

    if(!(n_ij == inner_segment_cv1.edge_int && n_ij == inner_segment_cv2.edge_int)){
        throw std::invalid_argument("get_close_scoop_sides: vertices are not on the same edge");
    }
    CombinatorialVertex aligned_cv1, aligned_cv2;
    align_cvs(inner_segment_cv1, inner_segment_cv2, aligned_cv1, aligned_cv2);

    auto kl = complement({i, j});
    int k = kl[0], l = kl[1];
    int smaller_ij_order = std::min(aligned_cv1.order, aligned_cv2.order);

    vector<int> close_sides;
    // check if it is a scoop on ijk side
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, k) &&
        smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, k) &&
        (smaller_ij_order - smeared_curve_nc.scc.get(i, j, k)) % 2 == 1){
        close_sides.push_back(k);
    }
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, l) &&
        smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, l) &&
        (smaller_ij_order - smeared_curve_nc.scc.get(i, j, l)) % 2 == 1){
        close_sides.push_back(l);
    }

    return close_sides;
}


// SPIRAL: walks the curve's own vertices and inserts an on-face scoop steiner
// in the middle of each edge segment. Never consults the NormalCoordinates.
CombFace
add_mid_scoop_vertices_spiral(
    const CombFace &spiral_comb_face
){
    vector<CombinatorialVertex> new_vertices;
    size_t num_vertices = spiral_comb_face.vertices.size();
    // add mid vertices on the open side of scoop segments
    for (size_t idx = 0; idx < num_vertices; idx++){
        CombinatorialVertex cv = spiral_comb_face.vertices[idx],
                            next_cv = spiral_comb_face.vertices[(idx+1) % num_vertices];
        new_vertices.push_back(cv);
        if (cv.same_edge_as(next_cv)){
            // this segment is along an edge; check if it is a scoop segment
            if(next_cv.k == -1) // this should hold construction from the boundary construction
                throw std::logic_error("add_mid_scoop_vertices: invalid segment edge with no k index");
            int other_vertex = next_cv.k;
            
            // min order of the aligned segment endpoints (same idiom as get_close_scoop_sides)
            CombinatorialVertex aligned_cv, aligned_next_cv;
            align_cvs(cv, next_cv, aligned_cv, aligned_next_cv);
            int min_order = std::min(aligned_cv.order, aligned_next_cv.order);
            
            CombinatorialVertex mid_cv(cv.i, cv.j, min_order, cv.edge_int, other_vertex);
            mid_cv.scoop_steiner = true;
            new_vertices.push_back(mid_cv);
        }
    }
    return CombFace(new_vertices, CombFaceType::SPIRAL);
}


// TUNNEL_QUAD: a single segment with 0, 1, or 2 close scoop sides.
static CombFace
add_mid_scoop_vertices_tunnel(
    const CombFace &tunnel_comb_face,
    const NormalCoordinates &smeared_curve_nc
){
    CombinatorialVertex cv1 = tunnel_comb_face.vertices[0],
                        cv2 = tunnel_comb_face.vertices[1];
    vector<int> close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2);

    if (close_scoop_sides.size() < 2){
        // no scoop on either side, no new polygons needed; stripe-to-stripe connection
        return tunnel_comb_face;
    }
    else if (close_scoop_sides.size() == 2){
        // Size 2 polygon basically
        int k = close_scoop_sides[0], l = close_scoop_sides[1];
        // add scoop steiner on both sides; no scoop_interior_steiner since there is no interior side
        CombinatorialVertex mid_on_face_cv1(cv1.i, cv1.j, cv1.order, cv1.edge_int, k);
        mid_on_face_cv1.scoop_steiner = true;
        CombinatorialVertex mid_on_face_cv2(cv1.i, cv1.j, cv1.order, cv1.edge_int, l);
        mid_on_face_cv2.scoop_steiner = true;
        return CombFace({cv1, mid_on_face_cv1, cv2, mid_on_face_cv2}, CombFaceType::TUNNEL_QUAD);
    }
    else {
        throw std::logic_error("add_mid_scoop_vertices: invalid number of close scoop sides");
    }
}


// Face-smeared types (SMEARED_QUAD, INTERIOR_HEXAGON, INTERIOR_CORNER_TYPE):
// walk the face's vertices; for each edge segment add either an on-face or a
// deep-interior scoop steiner depending on which side is the close scoop side.
static CombFace
add_mid_scoop_vertices_smeared(
    const CombFace &even_sum_comb_face,
    const NormalCoordinates &smeared_curve_nc
){
    size_t num_vertices = even_sum_comb_face.vertices.size();
    // the cv.k determines what face the segment is smeared on
    vector<CombinatorialVertex> new_cv_vertices;
    for (int idx = 0; idx < num_vertices; idx++){
        CombinatorialVertex cv1 = even_sum_comb_face.vertices[idx],
                            cv2 = even_sum_comb_face.vertices[(idx + 1) % num_vertices];
        new_cv_vertices.push_back(cv1);

        if (!cv1.same_edge_as(cv2)){
            // not an edge segment; no mid vertex needed
            continue;
        }
        vector<int> close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2);

        if (close_scoop_sides.size() == 2){
            throw std::logic_error("add_mid_scoop_vertices: invalid scenario should have been handled by the TUNNEL_QUAD case");
        }
        int i = cv1.i, j = cv1.j, k = cv1.k; // ijk is the smeared tet face
        // Add a mid vertex on the the closed side[s]

        // use min order so that the interpolator knows to just use order and order+1 for mid vertex construction
        // we already made sure cv1 and cv2 are on the same edge

        CombinatorialVertex aligned_cv1, aligned_cv2;
        align_cvs(cv1, cv2, aligned_cv1, aligned_cv2);
        int min_order = std::min({aligned_cv1.order, aligned_cv2.order});
        CombinatorialVertex mid_scoop_cv(aligned_cv1.i, aligned_cv1.j, min_order, aligned_cv1.edge_int, cv1.k);

        if (close_scoop_sides.size() == 0){
            // add mid vertex in the interior; it is simply pushed to interior of the tet a bit
            mid_scoop_cv.scoop_interior_steiner = true;
            new_cv_vertices.push_back(mid_scoop_cv);
        }
        else if (close_scoop_sides.size() == 1){
            int l = close_scoop_sides[0];
            // if l!=k, we could have any smear type, 
            // if l==k, then we have an interior corner type with multiple alternating scoops on the same edge
            if (l == i || l == j)
                throw std::logic_error("add_mid_scoop_vertices: invalid close scoop side");
            // use the scoop_steiner on the other face
            mid_scoop_cv.scoop_steiner = true;
            mid_scoop_cv.k = l; // change the k to the close scoop side
            new_cv_vertices.push_back(mid_scoop_cv);
        }
        else {
            throw std::logic_error("add_mid_scoop_vertices: invalid scenario with close scoop sides");
        }
    }

    return CombFace(new_cv_vertices, even_sum_comb_face.type);
}


CombFace
add_mid_scoop_vertices(
    const CombFace &even_sum_comb_face,
    const NormalCoordinates &smeared_curve_nc
){
    // different even-sum types are handled separately:
    // SPIRAL,
    // TUNNEL_QUAD,
    // SMEARED_QUAD / INTERIOR_HEXAGON / INTERIOR_CORNER_TYPE (face-smeared),
    // CORNER_TYPE_TRIANGLE (no mid vertices added)
    switch (even_sum_comb_face.type){
        case CombFaceType::SPIRAL:
            return add_mid_scoop_vertices_spiral(even_sum_comb_face);
        case CombFaceType::TUNNEL_QUAD:
            return add_mid_scoop_vertices_tunnel(even_sum_comb_face, smeared_curve_nc);
        case CombFaceType::SMEARED_QUAD:
        case CombFaceType::INTERIOR_HEXAGON:
        case CombFaceType::INTERIOR_CORNER_TYPE:
            return add_mid_scoop_vertices_smeared(even_sum_comb_face, smeared_curve_nc);
        default:
            throw std::logic_error("add_mid_scoop_vertices: unsupported comb face type");
    }
}
