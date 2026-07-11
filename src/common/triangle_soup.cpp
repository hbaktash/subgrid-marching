#include "common/triangle_soup.h"
#include "common/utils.h"


bool TriangleSoup::COMB_MERGE = false;


CombVertexKey
key_from_cv(const CombinatorialVertex &cv){
    CombVertexSigType type;
    int k_val = -1;
    if (cv.scoop_steiner){
        type = CombVertexSigType::SCOOP_FACE_STEINER;
        k_val = cv.k;
    } else if (cv.scoop_interior_steiner){
        type = CombVertexSigType::SCOOP_INTERIOR_STEINER;
    } else {
        type = CombVertexSigType::NORMAL;
    }
    int i = cv.i, j = cv.j, order = cv.order;
    if (i > j){
        std::swap(i, j);
        if (type == CombVertexSigType::NORMAL)
            order = cv.edge_int - order + 1;  // single intersection: flip order
        else
            order = cv.edge_int - order;       // steiner between [m,m+1]: flip to [n-m, n-m+1]
    }
    return CombVertexKey{i, j, order, cv.edge_int, k_val, type};
}


void
TriangleSoup::add_local_soup(const TriangleSoup &other){
    if (COMB_MERGE && !other.signatures.empty())
        add_local_soup_comb(other);
    else
        add_local_soup(other.faces, other.vertices);
}

void
TriangleSoup::add_local_soup(
    const vector<vector<size_t>> &local_faces,
    const vector<Vector3> &local_vertices
){
    size_t vertex_offset = vertices.size();
    for (const auto &v: local_vertices)
        vertices.push_back(v);
    for (const auto &face: local_faces){
        vector<size_t> offset_face;
        for (size_t idx: face)
            offset_face.push_back(idx + vertex_offset);
        faces.push_back(offset_face);
    }
}

void
TriangleSoup::add_local_soup_comb(const TriangleSoup &other){
    size_t n = other.vertices.size();
    vector<size_t> idx_map(n);
    for (size_t vi = 0; vi < n; vi++){
        const auto &key = other.signatures[vi];
        if (key.type == CombVertexSigType::STANDALONE){
            idx_map[vi] = vertices.size();
            vertices.push_back(other.vertices[vi]);
            signatures.push_back(key);
        } else {
            auto it = sig_to_idx.find(key);
            if (it != sig_to_idx.end()){
                idx_map[vi] = it->second;
            } else {
                idx_map[vi] = vertices.size();
                sig_to_idx[key] = idx_map[vi];
                vertices.push_back(other.vertices[vi]);
                signatures.push_back(key);
            }
        }
    }
    for (size_t fi = 0; fi < other.faces.size(); fi++){
        vector<size_t> mapped_face;
        for (size_t idx: other.faces[fi])
            mapped_face.push_back(idx_map[idx]);
        faces.push_back(mapped_face);
        if (!other.faces_per_edge.empty())
            faces_per_edge.push_back(other.faces_per_edge[fi]);
        if (!other.dual_positions.empty())
            dual_positions.push_back(other.dual_positions[fi]);
    }
}


void
TriangleSoup::make_vertex_standalone(int v){
    for (auto& sig : signatures){
        if (sig.i == v || sig.j == v){
            sig_to_idx.erase(sig);
            sig.type = CombVertexSigType::STANDALONE;
        }
    }
}


void
TriangleSoup::make_type_standalone(CombVertexSigType target_type){
    for (auto& sig : signatures){
        if (sig.type == target_type){
            sig_to_idx.erase(sig);
            sig.type = CombVertexSigType::STANDALONE;
        }
    }
}


void
TriangleSoup::remap_vertex_indices(const array<int, 4> &idx_map){
    if (signatures.empty()) return;
    for (auto &sig: signatures){
        if (sig.type == CombVertexSigType::STANDALONE) continue;
        sig.i = idx_map[sig.i];
        sig.j = idx_map[sig.j];
        if (sig.i > sig.j){
            std::swap(sig.i, sig.j);
            if (sig.type == CombVertexSigType::NORMAL)
                sig.order = sig.edge_int - sig.order + 1;
            else
                sig.order = sig.edge_int - sig.order;
        }
        if (sig.k >= 0) sig.k = idx_map[sig.k];
    }
    for (auto &fpe: faces_per_edge)
        for (auto &triplet: fpe)
            for (auto &v: triplet)
                v = (size_t)idx_map[v];
}


void
TriangleSoup::remap_orders(
    const array<vector<double>, 6> &edge_isect_ts,
    const EdgeOccupations &occ,
    bool count_occupied
){
    if (signatures.empty()) return;
    for (auto &sig: signatures){
        if (sig.type == CombVertexSigType::STANDALONE) continue;
        int edge_idx = edge_pair_to_index(sig.i, sig.j);
        int n = (int)edge_isect_ts[edge_idx].size();
        int count = 0;
        for (int k = 0; k < n; k++){
            bool occ_at_k = occ.get(sig.i, sig.j, k);
            bool should_count = count_occupied ? occ_at_k : !occ_at_k;
            if (should_count){
                count++;
                if (count == sig.order){
                    sig.order = k + 1;
                    sig.edge_int = n;
                    break;
                }
            }
        }
    }
}
