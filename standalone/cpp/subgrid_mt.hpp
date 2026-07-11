// subgrid_mt.hpp — self-contained, dependency-free port of the Subgrid Marching
// Tetrahedra per-tet reconstruction (primal + dual).
//
// STL only. Single header, copy-pastable. This is a FAITHFUL transcription of the
// repository core in src/subgrid_MT/ + src/nc/ + the pieces of common/utils it
// depends on. Do NOT change the algorithm here; if the core changes, re-port from
// the maintained core in src/subgrid_MT/.
//
// Two entry points, matching docs/construction_policy.md "Input":
//   subgrid_mt::subgrid_primal(tet_positions, tet_global_indices, edge_isect_ts)
//   subgrid_mt::subgrid_dual  (tet_positions, tet_global_indices, edge_isect_ts,
//                              edge_isect_normals, reg_alpha, project_duals)
//
// Both return a per-tet soup with global-remapped combinatorial signatures
// (already within-tet merged, as in the core). Global cross-tet merging is left to
// the caller; use signature_key() to dedup (see standalone/examples).

#pragma once

#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <initializer_list>
#include <cstddef>
#include <limits>
#include <iostream>
#include <cassert>

namespace subgrid_mt {

using std::vector;
using std::array;
using std::pair;
using std::tuple;

// ============================================================================
// Vector3  (minimal port of geometrycentral::Vector3)
// ============================================================================
struct Vec3 {
    double x = 0, y = 0, z = 0;

    static Vec3 zero() { return {0, 0, 0}; }
    static Vec3 undefined() {
        double n = std::numeric_limits<double>::quiet_NaN();
        return {n, n, n};
    }
    bool isDefined() const { return !(std::isnan(x) || std::isnan(y) || std::isnan(z)); }

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

    double norm2() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm2()); }
};

inline Vec3 operator*(double s, const Vec3& v) { return {v.x * s, v.y * s, v.z * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// ============================================================================
// tet constants + generic set/geometry helpers  (from common/utils)
// ============================================================================
inline const array<pair<int, int>, 6> ALL_TET_PAIRS = {
    pair<int, int>{0, 1}, pair<int, int>{0, 2}, pair<int, int>{0, 3},
    pair<int, int>{1, 2}, pair<int, int>{1, 3}, pair<int, int>{2, 3}};

inline const array<array<int, 3>, 4> ALL_TET_TRIPLETS = {
    array<int, 3>{1, 2, 3}, array<int, 3>{0, 2, 3},
    array<int, 3>{0, 1, 3}, array<int, 3>{0, 1, 2}};

inline size_t edge_pair_to_index(int i, int j) {
    if (i > j) std::swap(i, j);
    if (i == 0 && j == 1) return 0;
    if (i == 0 && j == 2) return 1;
    if (i == 0 && j == 3) return 2;
    if (i == 1 && j == 2) return 3;
    if (i == 1 && j == 3) return 4;
    if (i == 2 && j == 3) return 5;
    throw std::invalid_argument("edge_pair_to_index: Invalid edge pair");
}

inline vector<int> complement(std::initializer_list<int> inputs) {
    std::set<int> universe = {0, 1, 2, 3};
    for (int x : inputs) universe.erase(x);
    return vector<int>(universe.begin(), universe.end());
}
inline vector<int> complement(const vector<int>& inputs) {
    std::set<int> universe = {0, 1, 2, 3};
    for (int x : inputs) universe.erase(x);
    return vector<int>(universe.begin(), universe.end());
}

inline vector<pair<int, int>> all_pairs(vector<int> universe) {
    vector<pair<int, int>> pairs;
    for (size_t i = 0; i < universe.size(); i++)
        for (size_t j = i + 1; j < universe.size(); j++)
            pairs.push_back({universe[i], universe[j]});
    return pairs;
}
inline vector<pair<int, int>> all_pairs(array<int, 3> u) {
    return {{u[0], u[1]}, {u[0], u[2]}, {u[1], u[2]}};
}

inline bool triangle_inequality(int a, int b, int c, int& res_a, int& res_b, int& res_c) {
    int mx = std::max({a, b, c});
    bool feasible = (a + b + c) >= 2 * mx;
    res_a = 0; res_b = 0; res_c = 0;
    if (!feasible) {
        if (a == mx) { res_a = a - b - c; res_b = 0; res_c = 0; }
        else if (b == mx) { res_b = b - a - c; res_a = 0; res_c = 0; }
        else { res_c = c - a - b; res_a = 0; res_b = 0; }
    }
    return feasible;
}

inline bool segment_triangle_intersection(const Vec3& p0, const Vec3& p1, const Vec3& a,
                                          const Vec3& b, const Vec3& c, double& out_t) {
    const double EPS = 1e-12;
    Vec3 dir = p1 - p0;
    Vec3 edge1 = b - a;
    Vec3 edge2 = c - a;
    Vec3 pvec = cross(dir, edge2);
    double det = dot(edge1, pvec);
    if (std::abs(det) < EPS) return false;
    double invDet = 1.0 / det;
    Vec3 tvec = p0 - a;
    double u = dot(tvec, pvec) * invDet;
    if (u < -EPS || u > 1.0 + EPS) return false;
    Vec3 qvec = cross(tvec, edge1);
    double v = dot(dir, qvec) * invDet;
    if (v < -EPS || u + v > 1.0 + EPS) return false;
    double t_ray = dot(edge2, qvec) * invDet;
    double dir_len2 = dir.norm2();
    if (dir_len2 <= EPS) return false;
    double t_param = t_ray;
    if (t_param < -EPS || t_param > 1.0 + EPS) return false;
    out_t = std::max(0.0, std::min(1.0, t_param));
    return true;
}

inline bool segment_tet_intersection(const Vec3& p0, const Vec3& p1,
                                     const array<Vec3, 4>& tet_positions, double& out_t) {
    out_t = std::numeric_limits<double>::infinity();
    for (auto triplet : ALL_TET_TRIPLETS) {
        Vec3 a = tet_positions[triplet[0]];
        Vec3 b = tet_positions[triplet[1]];
        Vec3 c = tet_positions[triplet[2]];
        double t_face;
        if (segment_triangle_intersection(p0, p1, a, b, c, t_face))
            out_t = std::min(out_t, t_face);
    }
    return out_t != std::numeric_limits<double>::infinity();
}

inline bool is_inside_tet(const Vec3& p, const array<Vec3, 4>& tet_positions, double eps = 1e-12) {
    auto signed_tet_vol = [](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
        return dot(cross((b - a), (c - a)), (d - a));
    };
    double from_pd_vol = 0.0;
    for (auto triplet : ALL_TET_TRIPLETS)
        from_pd_vol += std::abs(signed_tet_vol(tet_positions[triplet[0]], tet_positions[triplet[1]],
                                               tet_positions[triplet[2]], p));
    double tet_vol = std::abs(signed_tet_vol(tet_positions[0], tet_positions[1],
                                             tet_positions[2], tet_positions[3]));
    return from_pd_vol <= tet_vol + eps;
}

inline vector<double> generate_uniform_isect_ts_single_edge(int edge_int, double t_max = 1.0) {
    vector<double> isect_ts;
    for (int i = 0; i < edge_int; i++)
        isect_ts.push_back(static_cast<double>(i + 1) * t_max / (double)(edge_int + 1));
    return isect_ts;
}

// ============================================================================
// Normal-coordinate containers  (from nc/normalCoord_containers.h)
// ============================================================================
struct EdgeInts {
    std::map<pair<int, int>, int> data;
    static pair<int, int> key(int i, int j) { if (i > j) std::swap(i, j); return {i, j}; }

    EdgeInts() = default;
    EdgeInts(array<int, 6> counts) {
        set(0, 1, counts[0]); set(0, 2, counts[1]); set(0, 3, counts[2]);
        set(1, 2, counts[3]); set(1, 3, counts[4]); set(2, 3, counts[5]);
    }
    EdgeInts(const array<vector<double>, 6>& isects) {
        set(0, 1, isects[0].size()); set(0, 2, isects[1].size()); set(0, 3, isects[2].size());
        set(1, 2, isects[3].size()); set(1, 3, isects[4].size()); set(2, 3, isects[5].size());
    }

    array<int, 6> to_array() const {
        return {get(0, 1), get(0, 2), get(0, 3), get(1, 2), get(1, 3), get(2, 3)};
    }
    int get(int i, int j) const {
        auto it = data.find(key(i, j));
        return it == data.end() ? 0 : it->second;
    }
    void set(int i, int j, size_t v) { data[key(i, j)] = (int)v; }
    int operator[](pair<int, int> p) const { return get(p.first, p.second); }
    int& operator[](pair<int, int> p) { return data[key(p.first, p.second)]; }
};

struct SemiCornerCuts {
    array<std::map<pair<int, int>, int>, 4> data;
    static pair<int, int> key(int j, int k) { if (j > k) std::swap(j, k); return {j, k}; }
    int get(int i, int j, int k) const {
        auto it = data[i].find(key(j, k));
        return it == data[i].end() ? 0 : it->second;
    }
    void set(int i, int j, int k, int v) { data[i][key(j, k)] = v; }
    const std::map<pair<int, int>, int>& operator[](size_t i) const { return data[i]; }
    std::map<pair<int, int>, int>& operator[](size_t i) { return data[i]; }
};

struct CornerCuts {
    int data[4] = {0, 0, 0, 0};
    int get(int i) const { return data[i]; }
    void set(int i, int v) { data[i] = v; }
    int operator[](int i) const { return data[i]; }
    int& operator[](int i) { return data[i]; }
};

struct DiagonalCuts {
    std::map<pair<int, int>, int> data;
    static pair<int, int> key(int i, int j) { if (i > j) std::swap(i, j); return {i, j}; }
    int get(int i, int j) const {
        auto it = data.find(key(i, j));
        return it == data.end() ? 0 : it->second;
    }
    void set(int i, int j, int v) {
        data[key(i, j)] = v;
        for (int k = 0; k < 4; k++) {
            if (k == i || k == j) continue;
            for (int l = k + 1; l < 4; l++) {
                if (l == i || l == j) continue;
                data[key(k, l)] = v;
            }
        }
    }
    int operator[](pair<int, int> p) const { return get(p.first, p.second); }
    int& operator[](pair<int, int> p) { return data[key(p.first, p.second)]; }
};

struct CombinatorialVertex {
    int i = -1, j = -1, order = -1, edge_int = 0;
    int k = -1;
    int l = -1;
    bool scoop_steiner = false;
    bool scoop_interior_steiner = false;

    CombinatorialVertex() = default;
    CombinatorialVertex(int i, int j, int order, int edge_int)
        : i(i), j(j), order(order), edge_int(edge_int) {}
    CombinatorialVertex(int i, int j, int order, int edge_int, int k)
        : i(i), j(j), order(order), edge_int(edge_int), k(k) {}
    CombinatorialVertex(int i, int j, int order, int edge_int, int k, int l)
        : i(i), j(j), order(order), edge_int(edge_int), k(k), l(l) {}

    bool operator==(const CombinatorialVertex& other) const {
        if ((i == other.i) && (j == other.j) && (order == other.order)) return true;
        if ((i == other.j) && (j == other.i) && (order == edge_int - other.order + 1)) return true;
        return false;
    }
    bool same_edge_as(const CombinatorialVertex& other) const {
        return ((i == other.i && j == other.j) || (i == other.j && j == other.i));
    }
};

struct DerivedVertex {
    vector<CombinatorialVertex> parents;
    vector<double> weights;

    DerivedVertex() = default;
    DerivedVertex(const vector<CombinatorialVertex>& parents)
        : parents(parents), weights(vector<double>(parents.size(), 1.0 / (double)parents.size())) {}
    DerivedVertex(const vector<CombinatorialVertex>& parents, const vector<double>& weights)
        : parents(parents), weights(weights) {}
    DerivedVertex(const DerivedVertex& dv1, const DerivedVertex& dv2, double t1) {
        parents = dv1.parents;
        parents.insert(parents.end(), dv2.parents.begin(), dv2.parents.end());
        auto weights1 = dv1.weights;
        for (double& w1 : weights1) w1 *= (1. - t1);
        weights = weights1;
        auto weights2 = dv2.weights;
        for (double& w2 : weights2) w2 *= t1;
        weights.insert(weights.end(), weights2.begin(), weights2.end());
    }
};

enum class CombFaceType {
    OPEN, CLOSED_NON_NORMAL, CLOSED_NORMAL,
    TRIANGLE, QUAD02, QUAD13,
    OCTAGON, SPIRAL,
    SMEARED_QUAD, INTERIOR_HEXAGON, INTERIOR_CORNER_TYPE,
    TUNNEL_QUAD
};

struct CombFace {
    vector<CombinatorialVertex> vertices;
    CombFaceType type;
    DerivedVertex derived_vertex;

    CombFace() = default;
    CombFace(const vector<CombinatorialVertex>& vertices, CombFaceType type)
        : vertices(vertices), type(type) {}
};

struct EdgeOccupations {
    std::map<pair<int, int>, vector<bool>> data;
    static pair<int, int> key(int i, int j) { if (i > j) std::swap(i, j); return {i, j}; }

    EdgeOccupations() = default;
    EdgeOccupations(const EdgeInts& edge_ints) {
        for (const auto& kv : edge_ints.data) data[kv.first] = vector<bool>(kv.second, false);
    }
    bool get(int i, int j, int index) const {
        if (i > j) index = data.at(key(i, j)).size() - 1 - index;
        return data.at(key(i, j))[index];
    }
    bool get(CombinatorialVertex cv) const {
        int i = cv.i, j = cv.j, index = cv.order - 1;
        if (i > j) index = data.at(key(i, j)).size() - 1 - index;
        return data.at(key(i, j))[index];
    }
    void set(int i, int j, int index, bool value) {
        if (i > j) index = data[key(i, j)].size() - 1 - index;
        data[key(i, j)][index] = value;
    }
    void set(CombinatorialVertex cv, bool value) {
        int i = cv.i, j = cv.j, index = cv.order - 1;
        if (i > j) index = data[key(i, j)].size() - 1 - index;
        data[key(i, j)][index] = value;
    }
    EdgeInts to_edge_ints(bool count_occupied) const {
        EdgeInts result;
        for (const auto& kv : data) {
            int count = 0;
            for (bool occupied : kv.second) {
                if (occupied && count_occupied) count++;
                if (!occupied && !count_occupied) count++;
            }
            result[kv.first] = count;
        }
        return result;
    }
};

// ============================================================================
// Combinatorial signature + TriangleSoup  (from common/triangle_soup.h/.cpp)
// ============================================================================
enum class CombVertexSigType { NORMAL, SCOOP_FACE_STEINER, SCOOP_INTERIOR_STEINER, STANDALONE };

struct CombVertexKey {
    int i = -1, j = -1;
    int order = -1;
    int edge_int = 0;  // NOT part of identity
    int k = -1;
    CombVertexSigType type = CombVertexSigType::STANDALONE;

    bool operator==(const CombVertexKey& o) const {
        return i == o.i && j == o.j && order == o.order && k == o.k && type == o.type;
    }
};

struct CombVertexKeyHash {
    size_t operator()(const CombVertexKey& k) const {
        size_t h = std::hash<int>{}(k.i);
        h ^= std::hash<int>{}(k.j) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(k.order) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(k.k) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(static_cast<int>(k.type)) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

inline CombVertexKey key_from_cv(const CombinatorialVertex& cv) {
    CombVertexSigType type;
    int k_val = -1;
    if (cv.scoop_steiner) { type = CombVertexSigType::SCOOP_FACE_STEINER; k_val = cv.k; }
    else if (cv.scoop_interior_steiner) { type = CombVertexSigType::SCOOP_INTERIOR_STEINER; }
    else { type = CombVertexSigType::NORMAL; }
    int i = cv.i, j = cv.j, order = cv.order;
    if (i > j) {
        std::swap(i, j);
        if (type == CombVertexSigType::NORMAL) order = cv.edge_int - order + 1;
        else order = cv.edge_int - order;
    }
    return CombVertexKey{i, j, order, cv.edge_int, k_val, type};
}

// Stable hashable string identity for cross-tet merging. Returns empty string for
// STANDALONE (which must NEVER merge — allocate a fresh vertex instead).
inline std::string signature_key(const CombVertexKey& s) {
    if (s.type == CombVertexSigType::STANDALONE) return std::string();
    return std::to_string(s.i) + "," + std::to_string(s.j) + "," + std::to_string(s.order) +
           "," + std::to_string(s.k) + "," + std::to_string(static_cast<int>(s.type));
}

struct TriangleSoup {
    // Set once before construction; read-only during. The port always runs true.
    static inline bool COMB_MERGE = true;

    vector<Vec3> vertices;
    vector<vector<size_t>> faces;
    vector<vector<array<size_t, 3>>> faces_per_edge;  // dual pipeline
    vector<Vec3> dual_positions;                       // dual pipeline
    std::unordered_map<CombVertexKey, size_t, CombVertexKeyHash> sig_to_idx;
    vector<CombVertexKey> signatures;

    void add_local_soup(const vector<vector<size_t>>& local_faces, const vector<Vec3>& local_vertices) {
        size_t vertex_offset = vertices.size();
        for (const auto& v : local_vertices) vertices.push_back(v);
        for (const auto& face : local_faces) {
            vector<size_t> offset_face;
            for (size_t idx : face) offset_face.push_back(idx + vertex_offset);
            faces.push_back(offset_face);
        }
    }

    void add_local_soup(const TriangleSoup& other) {
        if (COMB_MERGE && !other.signatures.empty()) add_local_soup_comb(other);
        else add_local_soup(other.faces, other.vertices);
    }

    void add_local_soup_comb(const TriangleSoup& other) {
        size_t n = other.vertices.size();
        vector<size_t> idx_map(n);
        for (size_t vi = 0; vi < n; vi++) {
            const auto& key = other.signatures[vi];
            if (key.type == CombVertexSigType::STANDALONE) {
                idx_map[vi] = vertices.size();
                vertices.push_back(other.vertices[vi]);
                signatures.push_back(key);
            } else {
                auto it = sig_to_idx.find(key);
                if (it != sig_to_idx.end()) {
                    idx_map[vi] = it->second;
                } else {
                    idx_map[vi] = vertices.size();
                    sig_to_idx[key] = idx_map[vi];
                    vertices.push_back(other.vertices[vi]);
                    signatures.push_back(key);
                }
            }
        }
        for (size_t fi = 0; fi < other.faces.size(); fi++) {
            vector<size_t> mapped_face;
            for (size_t idx : other.faces[fi]) mapped_face.push_back(idx_map[idx]);
            faces.push_back(mapped_face);
            if (!other.faces_per_edge.empty()) faces_per_edge.push_back(other.faces_per_edge[fi]);
            if (!other.dual_positions.empty()) dual_positions.push_back(other.dual_positions[fi]);
        }
    }

    void make_vertex_standalone(int v) {
        for (auto& sig : signatures) {
            if (sig.i == v || sig.j == v) {
                sig_to_idx.erase(sig);
                sig.type = CombVertexSigType::STANDALONE;
            }
        }
    }

    void make_type_standalone(CombVertexSigType target_type) {
        for (auto& sig : signatures) {
            if (sig.type == target_type) {
                sig_to_idx.erase(sig);
                sig.type = CombVertexSigType::STANDALONE;
            }
        }
    }

    void remap_vertex_indices(const array<int, 4>& idx_map) {
        if (signatures.empty()) return;
        for (auto& sig : signatures) {
            if (sig.type == CombVertexSigType::STANDALONE) continue;
            sig.i = idx_map[sig.i];
            sig.j = idx_map[sig.j];
            if (sig.i > sig.j) {
                std::swap(sig.i, sig.j);
                if (sig.type == CombVertexSigType::NORMAL) sig.order = sig.edge_int - sig.order + 1;
                else sig.order = sig.edge_int - sig.order;
            }
            if (sig.k >= 0) sig.k = idx_map[sig.k];
        }
        for (auto& fpe : faces_per_edge)
            for (auto& triplet : fpe)
                for (auto& v : triplet) v = (size_t)idx_map[v];
    }

    void remap_orders(const array<vector<double>, 6>& edge_isect_ts, const EdgeOccupations& occ,
                      bool count_occupied) {
        if (signatures.empty()) return;
        for (auto& sig : signatures) {
            if (sig.type == CombVertexSigType::STANDALONE) continue;
            int edge_idx = edge_pair_to_index(sig.i, sig.j);
            int n = (int)edge_isect_ts[edge_idx].size();
            int count = 0;
            for (int k = 0; k < n; k++) {
                bool occ_at_k = occ.get(sig.i, sig.j, k);
                bool should_count = count_occupied ? occ_at_k : !occ_at_k;
                if (should_count) {
                    count++;
                    if (count == sig.order) {
                        sig.order = k + 1;
                        sig.edge_int = n;
                        break;
                    }
                }
            }
        }
    }
};

// ============================================================================
// NormalCoordinates  (from nc/NC_solver.h/.cpp)
// ============================================================================
class NormalCoordinates {
public:
    EdgeInts edge_ints;
    SemiCornerCuts scc;
    CornerCuts cc;
    DiagonalCuts dc;
    SemiCornerCuts scoop;
    bool triangle_ineq = true, even_sum = true;

    NormalCoordinates() = default;
    NormalCoordinates(EdgeInts& edge_ints) : edge_ints(edge_ints) { compute_NC(); }
    NormalCoordinates(const EdgeInts& edge_ints) : edge_ints(edge_ints) { compute_NC(); }

    void compute_NC() {
        compute_semi_corner_cuts_and_scoops();
        if ((even_sum && triangle_ineq) || !even_sum) {
            SemiCornerCuts scc_reduced;
            compute_corner_cuts_and_reduced_semi_corner_cuts(scc_reduced);
            compute_diagonal_cuts(scc_reduced);
        } else if (even_sum) {
            EdgeInts reduced_edge_ints = scoop_reduction();
            NormalCoordinates reduced_NC;
            reduced_NC.edge_ints = reduced_edge_ints;
            reduced_NC.compute_NC();
            if (!reduced_NC.triangle_ineq || !reduced_NC.even_sum)
                throw std::logic_error("Error: compute_NC: scoop reduction did not yield feasible NC!\n");
            cc = reduced_NC.cc;
            dc = reduced_NC.dc;
        } else {
            return;
        }
    }

    void compute_semi_corner_cuts_and_scoops() {
        for (array<int, 3> triplet : ALL_TET_TRIPLETS) {
            int i = triplet[0], j = triplet[1], k = triplet[2];
            int res_ij, res_ik, res_jk;
            bool feasible = triangle_inequality(edge_ints[{i, j}], edge_ints[{i, k}], edge_ints[{j, k}],
                                                res_ij, res_ik, res_jk);
            bool ijk_even_sum = ((edge_ints[{i, j}] + edge_ints[{i, k}] + edge_ints[{j, k}]) % 2 == 0);
            triangle_ineq = feasible && triangle_ineq;
            even_sum = ijk_even_sum && even_sum;
            if (!ijk_even_sum && feasible) { res_ij = 1; res_ik = 1; res_jk = 1; }
            int M_ij = edge_ints[{i, j}] - res_ij, M_ik = edge_ints[{i, k}] - res_ik,
                M_jk = edge_ints[{j, k}] - res_jk;
            scc.set(i, j, k, std::max<int>((M_ij + M_ik - M_jk) / 2, 0));
            scc.set(j, i, k, std::max<int>((M_ij + M_jk - M_ik) / 2, 0));
            scc.set(k, i, j, std::max<int>((M_ik + M_jk - M_ij) / 2, 0));
            scoop.set(i, j, k, res_jk / 2);
            scoop.set(j, i, k, res_ik / 2);
            scoop.set(k, i, j, res_ij / 2);
        }
    }

    void compute_corner_cuts_and_reduced_semi_corner_cuts(SemiCornerCuts& scc_reduced) {
        cc[0] = std::min({scc.get(0, 1, 2), scc.get(0, 1, 3), scc.get(0, 2, 3)});
        cc[1] = std::min({scc.get(1, 0, 2), scc.get(1, 0, 3), scc.get(1, 2, 3)});
        cc[2] = std::min({scc.get(2, 0, 1), scc.get(2, 0, 3), scc.get(2, 1, 3)});
        cc[3] = std::min({scc.get(3, 0, 1), scc.get(3, 0, 2), scc.get(3, 1, 2)});
        for (int i = 0; i < 4; i++) {
            for (pair<int, int> jk_pair : all_pairs(complement({i}))) {
                int j = jk_pair.first, k = jk_pair.second;
                scc_reduced.set(i, j, k, scc.get(i, j, k) - cc[i]);
            }
        }
    }

    void compute_diagonal_cuts(const SemiCornerCuts& scc_reduced) {
        for (int j = 1; j < 4; j++) {
            vector<int> comp = complement({0, j});
            int k = comp[0], l = comp[1];
            dc.set(0, j, std::min({scc_reduced.get(k, 0, j), scc_reduced.get(l, 0, j),
                                   scc_reduced.get(j, k, l), scc_reduced.get(0, k, l)}));
        }
    }

    EdgeInts scoop_reduction() {
        EdgeOccupations edge_occupations(edge_ints);
        for (array<int, 3> triplet : ALL_TET_TRIPLETS) {
            int base, v1, v2;
            int i = triplet[0], j = triplet[1], k = triplet[2];
            if (scoop.get(i, j, k) > 0) { base = i; v1 = j; v2 = k; }
            else if (scoop.get(j, k, i) > 0) { base = j; v1 = k; v2 = i; }
            else if (scoop.get(k, i, j) > 0) { base = k; v1 = i; v2 = j; }
            else continue;
            int scoop_count = scoop.get(base, v1, v2);
            for (int s = 0; s < scoop_count; s++) {
                CombinatorialVertex cv_first(v1, v2, 2 * s + scc.get(v1, v2, base) + 1, edge_ints[{v1, v2}], base);
                CombinatorialVertex cv_second(v1, v2, 2 * s + scc.get(v1, v2, base) + 2, edge_ints[{v1, v2}], base);
                trace_forever(cv_first, edge_occupations);
                trace_forever(cv_second, edge_occupations);
            }
        }
        EdgeInts reduced_edge_ints;
        for (pair<int, int> ij_pair : all_pairs({0, 1, 2, 3})) {
            int i = ij_pair.first, j = ij_pair.second;
            int occ_count = 0;
            for (int idx = 0; idx < edge_ints[{i, j}]; idx++)
                if (!edge_occupations.get(i, j, idx)) occ_count++;
            reduced_edge_ints.set(i, j, occ_count);
        }
        return reduced_edge_ints;
    }

    CombinatorialVertex find_next_vertex(CombinatorialVertex v0) {
        if (even_sum == false)
            throw std::logic_error("combinatorial trace: even_sum is false, cannot trace scoops\n");
        int i = v0.i, j = v0.j;
        int k = complement({i, j, v0.k})[0];
        int order = v0.order;
        int n_ij = v0.edge_int;
        if (scc.get(i, j, k) + scc.get(j, i, k) + 2 * scoop.get(k, i, j) != n_ij)
            throw std::logic_error("combinatorial trace: inconsistent scoop and semicorner cut data\n");
        CombinatorialVertex v_next;
        if (order <= scc.get(i, j, k)) {
            v_next = CombinatorialVertex(i, k, order, edge_ints.get(i, k), j);
        } else if (order >= n_ij - scc.get(j, i, k) + 1) {
            v_next = CombinatorialVertex(j, k, n_ij - order + 1, edge_ints.get(j, k), i);
        } else {
            int next_order;
            if ((order - scc.get(i, j, k)) % 2 == 0) next_order = order - 1;
            else next_order = order + 1;
            v_next = CombinatorialVertex(i, j, next_order, n_ij, k);
        }
        return v_next;
    }

    void trace_forever(CombinatorialVertex v0, EdgeOccupations& edge_occupations) {
        CombinatorialVertex cv = v0;
        while (true) {
            if (edge_occupations.get(cv)) return;
            edge_occupations.set(cv, true);
            cv = find_next_vertex(cv);
        }
    }

    pair<int, int> get_d1_d2() const {
        int b1 = 0;
        int diag_types = 0;
        int d1 = 0, d2 = 0;
        for (int i = 1; i < 4; i++) {
            if (dc[{b1, i}] > 0) {
                diag_types++;
                if (d1 == 0) d1 = dc[{b1, i}];
                else d2 = dc[{b1, i}];
            }
        }
        if (diag_types > 2)
            throw std::logic_error("Error: get_d1_d2: more than 3 diagonal cut types, cannot determine d1, d2\n");
        return std::make_pair(d1, d2);
    }
};

// ============================================================================
// boundary curves  (from subgrid_MT/boundary_curve.cpp)
// ============================================================================
inline bool pick_residual_cv(int i, int j, int other, const NormalCoordinates& nc,
                             const array<size_t, 4>& global_inds, CombinatorialVertex& out_cv) {
    int nij = nc.edge_ints[{i, j}];
    int scoop_ij = nc.scoop.get(other, i, j);
    int scc_i = nc.scc.get(i, j, other);
    int scc_j = nc.scc.get(j, i, other);
    if (nij < scoop_ij * 2 + scc_i + scc_j) {
        throw std::logic_error("pick_residual_cv: inconsistent scoop and scc counts\n");
    } else if (nij == scoop_ij * 2 + scc_i + scc_j) {
        return false;
    } else if (nij == scoop_ij * 2 + scc_i + scc_j + 1) {
        if (scoop_ij == 0) {
            int order = nc.scc.get(i, j, other) + 1;
            out_cv = CombinatorialVertex(i, j, order, nij, other);
        } else {
            if (scoop_ij % 2 == 0) {
                int order = nc.scc.get(i, j, other) + scoop_ij + 1;
                out_cv = CombinatorialVertex(i, j, order, nij, other);
            } else {
                int order = nc.scc.get(i, j, other) + scoop_ij - 1 + 1;
                if (global_inds[i] > global_inds[j]) order += 2;
                out_cv = CombinatorialVertex(i, j, order, nij, other);
            }
        }
        return true;
    } else {
        throw std::logic_error("pick_residual_cv: more than 1 residual cv on this edge, cannot pick uniquely\n");
    }
}

inline CombinatorialVertex find_next_generic_vertex(const CombinatorialVertex& cv0,
                                                    const NormalCoordinates& nc,
                                                    const array<size_t, 4>& global_inds) {
    int i = cv0.i, j = cv0.j;
    int k = complement({i, j, cv0.k})[0];
    int order = cv0.order;
    int n_ij = cv0.edge_int;
    CombinatorialVertex v_next;
    if (order <= nc.scc.get(i, j, k)) {
        v_next = CombinatorialVertex(i, k, order, nc.edge_ints.get(i, k), j);
        return v_next;
    } else if (order >= n_ij - nc.scc.get(j, i, k) + 1) {
        v_next = CombinatorialVertex(j, k, n_ij - order + 1, nc.edge_ints.get(j, k), i);
        return v_next;
    }
    CombinatorialVertex cv_residual;
    bool has_residual = pick_residual_cv(i, j, k, nc, global_inds, cv_residual);
    int res_order = has_residual ? cv_residual.order : -1;
    if (order == res_order) {
        return cv_residual;
    } else {
        int next_order;
        if (order < res_order || !has_residual) {
            if ((order - nc.scc.get(i, j, k)) % 2 == 0) next_order = order - 1;
            else next_order = order + 1;
        } else {
            if ((order - nc.scc.get(i, j, k)) % 2 == 0) next_order = order + 1;
            else next_order = order - 1;
        }
        v_next = CombinatorialVertex(i, j, next_order, n_ij, k);
        return v_next;
    }
}

inline void trace_generic_curve(const CombinatorialVertex& start_cv, const NormalCoordinates& nc,
                                const array<size_t, 4>& global_inds, EdgeOccupations& occupations,
                                vector<CombinatorialVertex>& curve_vertices,
                                bool keep_occupation_track = true) {
    CombinatorialVertex current_cv = start_cv;
    if (keep_occupation_track) occupations.set(current_cv, true);
    curve_vertices.push_back(current_cv);
    while (true) {
        CombinatorialVertex next_cv = find_next_generic_vertex(current_cv, nc, global_inds);
        if (keep_occupation_track) occupations.set(next_cv, true);
        if (next_cv == current_cv || next_cv == start_cv) break;
        curve_vertices.push_back(next_cv);
        current_cv = next_cv;
    }
}

inline vector<CombinatorialVertex> pick_scoop_cvs(int i, int j, int other, const NormalCoordinates& nc,
                                                  const CombinatorialVertex& residual_cv) {
    assert((nc.edge_ints[{i, j}] >= 2 * nc.scoop.get(other, i, j) + nc.scc.get(i, j, other) + nc.scc.get(j, i, other)));
    int nij = nc.edge_ints[{i, j}];
    bool has_residual = residual_cv.i == -1;  // default cv
    int start_order = nc.scc.get(i, j, other) + 1, end_order = nij - nc.scc.get(j, i, other);
    vector<CombinatorialVertex> scoop_cvs;
    for (int order = start_order; order <= end_order; order++) {
        CombinatorialVertex cv(i, j, order, nij, other);
        if (!has_residual || !(cv == residual_cv)) scoop_cvs.push_back(cv);
    }
    return scoop_cvs;
}

inline void pick_residuals_and_scoop_cvs(const NormalCoordinates& NC, const array<size_t, 4>& tet_global_indices,
                                         vector<CombinatorialVertex>& residual_cvs,
                                         vector<CombinatorialVertex>& scoop_cvs) {
    for (const auto& ijk : ALL_TET_TRIPLETS) {
        for (auto ij : all_pairs(ijk)) {
            int i = ij.first, j = ij.second;
            int k = ijk[0] + ijk[1] + ijk[2] - i - j;
            CombinatorialVertex residual_cv;
            if (pick_residual_cv(i, j, k, NC, tet_global_indices, residual_cv))
                residual_cvs.push_back(residual_cv);
            auto ij_scoop_cvs = pick_scoop_cvs(i, j, k, NC, residual_cv);
            scoop_cvs.insert(scoop_cvs.end(), ij_scoop_cvs.begin(), ij_scoop_cvs.end());
        }
    }
}

inline tuple<vector<CombFace>, vector<CombFace>, vector<CombFace>>
boundary_comb_curves(const array<size_t, 4> tet_global_indices, const array<int, 6> edge_isect_counts,
                     EdgeOccupations& out_edge_occupations, bool populate_normal_curves) {
    EdgeInts edge_ints{edge_isect_counts};
    NormalCoordinates NC;
    NC.edge_ints = edge_ints;
    NC.compute_semi_corner_cuts_and_scoops();

    vector<CombinatorialVertex> residual_cvs, scoop_cvs;
    pick_residuals_and_scoop_cvs(NC, tet_global_indices, residual_cvs, scoop_cvs);

    EdgeOccupations edge_occupations(edge_ints);
    vector<CombFace> open_curves;
    for (const auto& res_cv : residual_cvs) {
        if (!edge_occupations.get(res_cv)) {
            vector<CombinatorialVertex> curve_vertices;
            trace_generic_curve(res_cv, NC, tet_global_indices, edge_occupations, curve_vertices);
            open_curves.emplace_back(curve_vertices, CombFaceType::OPEN);
        }
    }

    vector<CombFace> scoop_curves;
    for (const auto& scoop_cv : scoop_cvs) {
        if (!edge_occupations.get(scoop_cv)) {
            vector<CombinatorialVertex> curve_vertices;
            trace_generic_curve(scoop_cv, NC, tet_global_indices, edge_occupations, curve_vertices);
            scoop_curves.emplace_back(curve_vertices, CombFaceType::CLOSED_NON_NORMAL);
        }
    }

    out_edge_occupations = edge_occupations;

    vector<CombFace> normal_curves;
    if (populate_normal_curves) {
        for (auto ij : ALL_TET_PAIRS) {
            int i = ij.first, j = ij.second;
            for (int order = 1; order <= edge_ints[{i, j}]; order++) {
                CombinatorialVertex cv(i, j, order, edge_ints[{i, j}]);
                if (!edge_occupations.get(cv)) {
                    vector<CombinatorialVertex> curve_vertices;
                    trace_generic_curve(cv, NC, tet_global_indices, edge_occupations, curve_vertices);
                    normal_curves.emplace_back(curve_vertices, CombFaceType::CLOSED_NORMAL);
                }
            }
        }
    }
    return {open_curves, scoop_curves, normal_curves};
}

inline array<size_t, 3> segment_tet_face(const CombinatorialVertex& v_a, const CombinatorialVertex& v_b) {
    if (v_a.same_edge_as(v_b)) return {(size_t)v_b.i, (size_t)v_b.j, (size_t)v_b.k};
    int third = (v_a.i != v_b.i && v_a.i != v_b.j) ? v_a.i : v_a.j;
    return {(size_t)v_b.i, (size_t)v_b.j, (size_t)third};
}

inline vector<array<size_t, 3>> curve_segment_faces(const CombFace& comb_face) {
    const auto& verts = comb_face.vertices;
    bool is_closed = (comb_face.type != CombFaceType::OPEN);
    size_t n = verts.size();
    size_t n_segments = is_closed ? n : n - 1;
    vector<array<size_t, 3>> faces;
    faces.reserve(n_segments);
    for (size_t i = 0; i < n_segments; ++i) faces.push_back(segment_tet_face(verts[i], verts[(i + 1) % n]));
    return faces;
}

// ============================================================================
// cv interpolation  (from subgrid_MT/cv_interpolation.cpp)
// ============================================================================
inline Vec3 interpolate_comb_vertex(const CombinatorialVertex& cv, const array<Vec3, 4>& tet_positions,
                                    const array<vector<double>, 6>& edge_isect_ts, double bulge) {
    int i = cv.i, j = cv.j, order = cv.order, n_ij = cv.edge_int;
    int edge_idx = edge_pair_to_index(i, j);
    assert((order <= (int)edge_isect_ts[edge_idx].size()) && (cv.edge_int == (int)edge_isect_ts[edge_idx].size()));
    double t;
    if (i < j) {
        t = edge_isect_ts[edge_idx][order - 1];
    } else {
        int tmp_order = n_ij - order + 1;
        t = 1 - edge_isect_ts[edge_idx][tmp_order - 1];
    }
    Vec3 pos = (1 - t) * tet_positions[i] + t * tet_positions[j];
    if (!cv.scoop_interior_steiner && !cv.scoop_steiner) {
        return pos;
    } else {
        if (cv.scoop_interior_steiner && cv.scoop_steiner)
            throw std::invalid_argument("interpolate_comb_vertex: vertex cannot be both scoop_steiner and scoop_interior_steiner");
        double next_t;
        int next_order = order + 1;
        if (i < j) {
            next_t = edge_isect_ts[edge_idx][next_order - 1];
        } else {
            int tmp_next_order = n_ij - next_order + 1;
            next_t = 1 - edge_isect_ts[edge_idx][tmp_next_order - 1];
        }
        Vec3 next_pos = (1 - next_t) * tet_positions[i] + next_t * tet_positions[j];
        Vec3 mid_pos = 0.5 * (next_pos + pos);
        if (cv.scoop_interior_steiner) {
            Vec3 mid_other_edge = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3] - tet_positions[i] - tet_positions[j]) / 2.0;
            return (1 - bulge) * mid_pos + bulge * mid_other_edge;
        } else {
            int k = cv.k;
            Vec3 face_interior_point = (2.0 * tet_positions[k] + tet_positions[i] + tet_positions[j]) / 4.0;
            return (1 - bulge) * mid_pos + bulge * face_interior_point;
        }
    }
}

inline Vec3 interpolate_derived_vertex(const DerivedVertex& dv, const array<Vec3, 4>& tet_positions,
                                       const array<vector<double>, 6>& edge_isect_ts) {
    Vec3 pos = Vec3::zero();
    for (size_t p = 0; p < dv.parents.size(); p++) {
        Vec3 parent_pos = interpolate_comb_vertex(dv.parents[p], tet_positions, edge_isect_ts, -1.0);
        pos += dv.weights[p] * parent_pos;
    }
    return pos;
}

namespace detail {

inline TriangleSoup embed_triangle_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                        const array<vector<double>, 6>& edge_isect_ts) {
    if (comb_face.vertices.size() != 3)
        throw std::invalid_argument("embed_triangle_face: CombFace of type TRIANGLE must have exactly 3 vertices");
    TriangleSoup soup;
    for (const auto& cv : comb_face.vertices) {
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {{0, 1, 2}};
    return soup;
}

inline TriangleSoup embed_quad_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                    const array<vector<double>, 6>& edge_isect_ts) {
    if (comb_face.vertices.size() != 4)
        throw std::invalid_argument("embed_quad_face: CombFace of type QUAD02 or QUAD13 must have exactly 4 vertices");
    vector<size_t> tri1 = comb_face.type == CombFaceType::QUAD02 ? vector<size_t>{0, 1, 2} : vector<size_t>{0, 1, 3};
    vector<size_t> tri2 = comb_face.type == CombFaceType::QUAD02 ? vector<size_t>{2, 3, 0} : vector<size_t>{1, 2, 3};
    TriangleSoup soup;
    for (const auto& cv : comb_face.vertices) {
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {tri1, tri2};
    return soup;
}

inline TriangleSoup embed_octagon_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                       const array<vector<double>, 6>& edge_isect_ts) {
    if (comb_face.vertices.size() != 8)
        throw std::invalid_argument("embed_octagon_face: CombFace of type OCTAGON must have exactly 8 vertices");
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    for (size_t i = 0; i < n; i++) {
        soup.vertices.push_back(interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, -1.0));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(comb_face.vertices[i]));
        soup.faces.push_back({n, i, (i + 1) % n});
    }
    soup.vertices.push_back(interpolate_derived_vertex(comb_face.derived_vertex, tet_positions, edge_isect_ts));
    if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
    return soup;
}

inline TriangleSoup embed_spiral_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                      const array<vector<double>, 6>& edge_isect_ts, double bulge) {
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    Vec3 mean_pos = Vec3::zero();
    for (size_t i = 0; i < n; i++) {
        Vec3 pos = interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, bulge);
        soup.vertices.push_back(pos);
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(comb_face.vertices[i]));
        mean_pos += pos;
        soup.faces.push_back({n, i, (i + 1) % n});
    }
    mean_pos /= (double)n;
    soup.vertices.push_back(mean_pos);
    if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
    return soup;
}

inline TriangleSoup embed_interior_hexagon_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                                const array<vector<double>, 6>& edge_isect_ts, double bulge) {
    if (comb_face.vertices.size() != 9)
        throw std::invalid_argument("embed_interior_hexagon_face: CombFace of type INTERIOR_HEXAGON with mid_scoop_vertices=true must have exactly 9 vertices (6 original + 3 mid scoop)");
    TriangleSoup soup;
    for (const auto& cv : comb_face.vertices) {
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    soup.faces = {{8, 2, 5}, {2, 8, 0}, {2, 0, 1}, {5, 2, 3}, {5, 3, 4}, {8, 5, 6}, {8, 6, 7}};
    return soup;
}

inline Vec3 interior_corner_type_bulge_steiner(const CombFace& comb_face, const vector<Vec3>& positions,
                                               const array<Vec3, 4>& tet_positions, double bulge) {
    CombinatorialVertex cv0 = comb_face.vertices[0], cv2 = comb_face.vertices[2];
    int common_i = (cv0.i == cv2.i || cv0.i == cv2.j) ? cv0.i : (cv0.j == cv2.i || cv0.j == cv2.j) ? cv0.j : -1;
    if (common_i == -1)
        throw std::invalid_argument("embed_corner_type_triangle_face: first two vertices must be on the same edge");
    int l = complement({cv0.i, cv0.j, cv2.i, cv2.j})[0];
    Vec3 mid_pos = 0.5 * (positions[0] + positions[2]);
    Vec3 mid_other_edge = (tet_positions[l] + tet_positions[common_i]) / 2.0;
    return (1 - bulge) * mid_pos + bulge * mid_other_edge;
}

inline TriangleSoup embed_interior_corner_type_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                                    const array<vector<double>, 6>& edge_isect_ts,
                                                    bool mid_scoop_vertices, double bulge) {
    TriangleSoup soup;
    for (const auto& cv : comb_face.vertices) {
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    size_t n = soup.vertices.size();
    if (!mid_scoop_vertices) {
        for (size_t i = 1; i < n - 1; i++) soup.faces.push_back({0, i, i + 1});
    } else {
        Vec3 bulge_steiner = interior_corner_type_bulge_steiner(comb_face, soup.vertices, tet_positions, bulge);
        soup.vertices.push_back(bulge_steiner);
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(CombVertexKey{-1, -1, -1, 0, -1, CombVertexSigType::STANDALONE});
        for (size_t i = 0; i < n; i++) soup.faces.push_back({n, i, (i + 1) % n});
    }
    return soup;
}

inline TriangleSoup embed_smeared_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                       const array<vector<double>, 6>& edge_isect_ts, double bulge) {
    size_t n = comb_face.vertices.size();
    TriangleSoup soup;
    for (const auto& cv : comb_face.vertices) {
        soup.vertices.push_back(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge));
        if (TriangleSoup::COMB_MERGE) soup.signatures.push_back(key_from_cv(cv));
    }
    for (size_t i = 2; i < n; i++) soup.faces.push_back({1, i, (i + 1) % n});
    return soup;
}

}  // namespace detail

inline TriangleSoup embed_comb_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                    const array<vector<double>, 6>& edge_isect_ts, bool mid_scoop_vertices,
                                    double bulge) {
    if (comb_face.vertices.size() < 3)
        throw std::invalid_argument("embed_comb_face: CombFace expects at least 3 vertices. Tunnels should have already received the extra cv's if needed.");
    if (comb_face.type == CombFaceType::TRIANGLE)
        return detail::embed_triangle_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::QUAD02 || comb_face.type == CombFaceType::QUAD13)
        return detail::embed_quad_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::OCTAGON)
        return detail::embed_octagon_face(comb_face, tet_positions, edge_isect_ts);
    else if (comb_face.type == CombFaceType::SPIRAL)
        return detail::embed_spiral_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else if (comb_face.type == CombFaceType::INTERIOR_CORNER_TYPE)
        return detail::embed_interior_corner_type_face(comb_face, tet_positions, edge_isect_ts, mid_scoop_vertices, bulge);
    else if (mid_scoop_vertices && comb_face.type == CombFaceType::INTERIOR_HEXAGON)
        return detail::embed_interior_hexagon_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else if (comb_face.type == CombFaceType::SMEARED_QUAD || comb_face.type == CombFaceType::INTERIOR_HEXAGON ||
             comb_face.type == CombFaceType::TUNNEL_QUAD)
        return detail::embed_smeared_face(comb_face, tet_positions, edge_isect_ts, bulge);
    else
        throw std::invalid_argument("embed_comb_face: CombFace type not supported");
}

// ============================================================================
// simplicial embedding  (from subgrid_MT/simplicial_embedding.cpp)
// ============================================================================
inline void align_cvs(const CombinatorialVertex& cv1, const CombinatorialVertex& cv2,
                      CombinatorialVertex& aligned_cv1, CombinatorialVertex& aligned_cv2) {
    if (!cv1.same_edge_as(cv2)) throw std::invalid_argument("align_cvs: vertices are not on the same edge");
    if (cv1.edge_int != cv2.edge_int)
        throw std::invalid_argument("align_cvs: vertices do not have the same edge intersection count");
    if (cv1.i == cv2.i && cv1.j == cv2.j) {
        aligned_cv1 = cv1;
        aligned_cv2 = cv2;
    } else if (cv1.i == cv2.j && cv1.j == cv2.i) {
        int n_ij = cv1.edge_int;
        aligned_cv1 = CombinatorialVertex(cv1.i, cv1.j, cv1.order, n_ij, cv1.k);
        aligned_cv2 = CombinatorialVertex(cv1.i, cv1.j, n_ij - cv2.order + 1, n_ij, cv2.k);
    } else {
        throw std::logic_error("align_cvs: invalid vertices");
    }
}

inline vector<int> get_close_scoop_sides(const NormalCoordinates& smeared_curve_nc,
                                         const CombinatorialVertex& inner_segment_cv1,
                                         const CombinatorialVertex& inner_segment_cv2) {
    if (!inner_segment_cv1.same_edge_as(inner_segment_cv2))
        throw std::invalid_argument("get_close_scoop_sides: vertices are not on the same edge");
    int i = inner_segment_cv1.i, j = inner_segment_cv1.j, n_ij = smeared_curve_nc.edge_ints[{i, j}];
    if (!(n_ij == inner_segment_cv1.edge_int && n_ij == inner_segment_cv2.edge_int))
        throw std::invalid_argument("get_close_scoop_sides: vertices are not on the same edge");
    CombinatorialVertex aligned_cv1, aligned_cv2;
    align_cvs(inner_segment_cv1, inner_segment_cv2, aligned_cv1, aligned_cv2);
    auto kl = complement({i, j});
    int k = kl[0], l = kl[1];
    int smaller_ij_order = std::min(aligned_cv1.order, aligned_cv2.order);
    vector<int> close_sides;
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, k) &&
        smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, k) &&
        (smaller_ij_order - smeared_curve_nc.scc.get(i, j, k)) % 2 == 1)
        close_sides.push_back(k);
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, l) &&
        smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, l) &&
        (smaller_ij_order - smeared_curve_nc.scc.get(i, j, l)) % 2 == 1)
        close_sides.push_back(l);
    return close_sides;
}

inline CombFace add_mid_scoop_vertices_spiral(const CombFace& spiral_comb_face) {
    vector<CombinatorialVertex> new_vertices;
    size_t num_vertices = spiral_comb_face.vertices.size();
    for (size_t idx = 0; idx < num_vertices; idx++) {
        CombinatorialVertex cv = spiral_comb_face.vertices[idx],
                            next_cv = spiral_comb_face.vertices[(idx + 1) % num_vertices];
        new_vertices.push_back(cv);
        if (cv.same_edge_as(next_cv)) {
            if (next_cv.k == -1)
                throw std::logic_error("add_mid_scoop_vertices: invalid segment edge with no k index");
            int other_vertex = next_cv.k;
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

namespace detail {

inline CombFace add_mid_scoop_vertices_tunnel(const CombFace& tunnel_comb_face,
                                              const NormalCoordinates& smeared_curve_nc) {
    CombinatorialVertex cv1 = tunnel_comb_face.vertices[0], cv2 = tunnel_comb_face.vertices[1];
    vector<int> close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2);
    if (close_scoop_sides.size() < 2) {
        return tunnel_comb_face;
    } else if (close_scoop_sides.size() == 2) {
        int k = close_scoop_sides[0], l = close_scoop_sides[1];
        CombinatorialVertex mid_on_face_cv1(cv1.i, cv1.j, cv1.order, cv1.edge_int, k);
        mid_on_face_cv1.scoop_steiner = true;
        CombinatorialVertex mid_on_face_cv2(cv1.i, cv1.j, cv1.order, cv1.edge_int, l);
        mid_on_face_cv2.scoop_steiner = true;
        return CombFace({cv1, mid_on_face_cv1, cv2, mid_on_face_cv2}, CombFaceType::TUNNEL_QUAD);
    } else {
        throw std::logic_error("add_mid_scoop_vertices: invalid number of close scoop sides");
    }
}

inline CombFace add_mid_scoop_vertices_smeared(const CombFace& even_sum_comb_face,
                                               const NormalCoordinates& smeared_curve_nc) {
    size_t num_vertices = even_sum_comb_face.vertices.size();
    vector<CombinatorialVertex> new_cv_vertices;
    for (int idx = 0; idx < (int)num_vertices; idx++) {
        CombinatorialVertex cv1 = even_sum_comb_face.vertices[idx],
                            cv2 = even_sum_comb_face.vertices[(idx + 1) % num_vertices];
        new_cv_vertices.push_back(cv1);
        if (!cv1.same_edge_as(cv2)) continue;
        vector<int> close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2);
        if (close_scoop_sides.size() == 2)
            throw std::logic_error("add_mid_scoop_vertices: invalid scenario should have been handled by the TUNNEL_QUAD case");
        int i = cv1.i, j = cv1.j;
        CombinatorialVertex aligned_cv1, aligned_cv2;
        align_cvs(cv1, cv2, aligned_cv1, aligned_cv2);
        int min_order = std::min({aligned_cv1.order, aligned_cv2.order});
        CombinatorialVertex mid_scoop_cv(aligned_cv1.i, aligned_cv1.j, min_order, aligned_cv1.edge_int, cv1.k);
        if (close_scoop_sides.size() == 0) {
            mid_scoop_cv.scoop_interior_steiner = true;
            new_cv_vertices.push_back(mid_scoop_cv);
        } else if (close_scoop_sides.size() == 1) {
            int l = close_scoop_sides[0];
            if (l == i || l == j) throw std::logic_error("add_mid_scoop_vertices: invalid close scoop side");
            mid_scoop_cv.scoop_steiner = true;
            mid_scoop_cv.k = l;
            new_cv_vertices.push_back(mid_scoop_cv);
        } else {
            throw std::logic_error("add_mid_scoop_vertices: invalid scenario with close scoop sides");
        }
    }
    return CombFace(new_cv_vertices, even_sum_comb_face.type);
}

}  // namespace detail

inline CombFace add_mid_scoop_vertices(const CombFace& even_sum_comb_face,
                                       const NormalCoordinates& smeared_curve_nc) {
    switch (even_sum_comb_face.type) {
        case CombFaceType::SPIRAL:
            return add_mid_scoop_vertices_spiral(even_sum_comb_face);
        case CombFaceType::TUNNEL_QUAD:
            return detail::add_mid_scoop_vertices_tunnel(even_sum_comb_face, smeared_curve_nc);
        case CombFaceType::SMEARED_QUAD:
        case CombFaceType::INTERIOR_HEXAGON:
        case CombFaceType::INTERIOR_CORNER_TYPE:
            return detail::add_mid_scoop_vertices_smeared(even_sum_comb_face, smeared_curve_nc);
        default:
            throw std::logic_error("add_mid_scoop_vertices: unsupported comb face type");
    }
}

// ============================================================================
// evensum construction  (from subgrid_MT/evensum_construction.cpp)
// ============================================================================
namespace detail {

inline int get_interior_vertex(const vector<int>& separated_vs) {
    int interior_vertex;
    if (separated_vs.size() == 4 || separated_vs.size() == 0) interior_vertex = -1;
    else if (separated_vs.size() == 3) { vector<int> complement_vs = complement(separated_vs); interior_vertex = complement_vs[0]; }
    else if (separated_vs.size() == 1) interior_vertex = separated_vs[0];
    else throw std::logic_error("get_interior_vertex: Invalid number of separated vertices");
    if (interior_vertex < -1 || interior_vertex >= 4)
        throw std::logic_error("get_interior_vertex: computed interior vertex out of range [-1, 3]");
    return interior_vertex;
}

inline pair<array<vector<double>, 6>, EdgeOccupations>
extract_curve_edge_isect_ts(const vector<CombinatorialVertex>& curve_cvs,
                            const array<vector<double>, 6>& global_edge_isect_ts) {
    EdgeInts curve_edge_ints(global_edge_isect_ts);
    EdgeOccupations curve_edge_occupations(curve_edge_ints);
    for (const auto& cv : curve_cvs) curve_edge_occupations.set(cv, true);
    array<vector<double>, 6> curve_edge_isect_ts;
    for (int e = 0; e < 6; e++) {
        int i = ALL_TET_PAIRS[e].first, j = ALL_TET_PAIRS[e].second;
        for (size_t k = 0; k < global_edge_isect_ts[e].size(); k++)
            if (curve_edge_occupations.get(i, j, k)) curve_edge_isect_ts[e].push_back(global_edge_isect_ts[e][k]);
    }
    return {curve_edge_isect_ts, curve_edge_occupations};
}

inline TriangleSoup emit_comb_face(const CombFace& comb_face, const array<Vec3, 4>& tet_positions,
                                   const array<vector<double>, 6>& edge_isect_ts, bool scoop_mid_vertices,
                                   double scoop_mid_vertex_bulge) {
    return embed_comb_face(comb_face, tet_positions, edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge);
}

}  // namespace detail

inline EdgeInts combface_to_edgeints(const CombFace& comb_face) {
    EdgeInts result;
    for (auto& cv : comb_face.vertices) { int i = cv.i, j = cv.j; result[{i, j}]++; }
    return result;
}

inline vector<int> separated_vertices(const EdgeInts& cruve_edge_ints) {
    vector<int> result = {0};
    for (int i = 1; i < 4; i++)
        if (cruve_edge_ints[{0, i}] % 2 == 0) result.push_back(i);
    return result;
}

inline vector<CombFace> smeared_quads_construction(const NormalCoordinates& smeared_curve_nc,
                                                   const int& interior_vertex) {
    vector<CombFace> faces;
    for (const auto& ijk : ALL_TET_TRIPLETS) {
        for (int i : ijk) {
            int j = ijk[0] == i ? ijk[1] : ijk[0];
            int k = ijk[0] + ijk[1] + ijk[2] - i - j;
            int n_ij = smeared_curve_nc.edge_ints[{i, j}], n_ik = smeared_curve_nc.edge_ints[{i, k}];
            bool i_is_interior = (interior_vertex == i);
            for (int segment_order = 1; segment_order <= smeared_curve_nc.scc.get(i, j, k); segment_order++) {
                bool construct = (i_is_interior && segment_order % 2 == 1) || (!i_is_interior && segment_order % 2 == 0);
                if (construct && segment_order > 1) {
                    CombinatorialVertex cvij_1(i, j, segment_order - 1, n_ij, k);
                    CombinatorialVertex cvij_2(i, j, segment_order, n_ij, k);
                    CombinatorialVertex cvik_2(i, k, segment_order, n_ik, j);
                    CombinatorialVertex cvik_1(i, k, segment_order - 1, n_ik, j);
                    faces.push_back(CombFace({cvij_1, cvij_2, cvik_2, cvik_1}, CombFaceType::SMEARED_QUAD));
                }
            }
        }
    }
    return faces;
}

inline vector<CombFace> tunnel_construction(const NormalCoordinates& smeared_curve_nc, const int& interior_vertex) {
    vector<CombFace> faces;
    for (const auto& ij : ALL_TET_PAIRS) {
        int i = ij.first, j = ij.second;
        int n_ij = smeared_curve_nc.edge_ints[{i, j}];
        for (int segment_order = 2; segment_order <= n_ij; segment_order++) {
            bool i_is_interior = (interior_vertex == i);
            bool construct = (i_is_interior && segment_order % 2 == 1) || (!i_is_interior && segment_order % 2 == 0);
            if (construct) {
                CombinatorialVertex cvij_1(i, j, segment_order - 1, n_ij);
                CombinatorialVertex cvij_2(i, j, segment_order, n_ij);
                faces.push_back(CombFace({cvij_1, cvij_2}, CombFaceType::TUNNEL_QUAD));
            }
        }
    }
    return faces;
}

inline vector<CombFace> interior_region_construction(const NormalCoordinates& smeared_curve_nc,
                                                     const int& interior_vertex) {
    if (interior_vertex == -1) {
        for (auto& ijk : ALL_TET_TRIPLETS) {
            int i = ijk[0], j = ijk[1], k = ijk[2];
            if (smeared_curve_nc.scc.get(i, j, k) % 2 == 1) {
                if (smeared_curve_nc.scc.get(i, j, k) != 1 || smeared_curve_nc.scc.get(j, i, k) != 1 ||
                    smeared_curve_nc.scc.get(k, i, j) != 1)
                    throw std::logic_error("interior_region_construction: invalid contractible curve; should have exactly one segment on each face");
                CombinatorialVertex cv_ij(i, j, 1, 2, k);
                CombinatorialVertex cv_ik(i, k, 1, 2, j);
                CombinatorialVertex cv_ji(j, i, 1, 2, k);
                CombinatorialVertex cv_jk(j, k, 1, 2, i);
                CombinatorialVertex cv_ki(k, i, 1, 2, j);
                CombinatorialVertex cv_kj(k, j, 1, 2, i);
                return {CombFace({cv_ij, cv_ik, cv_ki, cv_kj, cv_jk, cv_ji}, CombFaceType::INTERIOR_HEXAGON)};
            }
        }
        return {};
    } else {
        int i = interior_vertex;
        int j, k;
        int l;
        auto jkl = complement({i});
        if (smeared_curve_nc.scc.get(i, jkl[0], jkl[1]) == 0) { j = jkl[0]; k = jkl[1]; }
        else if (smeared_curve_nc.scc.get(i, jkl[0], jkl[2]) == 0) { j = jkl[0]; k = jkl[2]; }
        else if (smeared_curve_nc.scc.get(i, jkl[1], jkl[2]) == 0) { j = jkl[1]; k = jkl[2]; }
        else throw std::invalid_argument("interior_region_construction: invalid corner curve; should have at least one open edge");
        l = 6 - i - j - k;
        CombinatorialVertex cv_ij(i, j, 1, smeared_curve_nc.edge_ints[{i, j}], k);
        CombinatorialVertex cv_il(i, l, 1, smeared_curve_nc.edge_ints[{i, l}], j);
        CombinatorialVertex cv_ik(i, k, 1, smeared_curve_nc.edge_ints[{i, k}], j);
        CombinatorialVertex cv_kj(k, j, smeared_curve_nc.scc.get(k, j, i), smeared_curve_nc.edge_ints[{k, j}], i);
        vector<CombinatorialVertex> interior_region_vertices = {cv_ij, cv_il, cv_ik, cv_kj};
        int n_kj = smeared_curve_nc.edge_ints[{k, j}];
        for (int scoop_idx = 0; scoop_idx < smeared_curve_nc.scoop.get(i, j, k); scoop_idx++) {
            CombinatorialVertex scoop_cv1(k, j, 2 * scoop_idx + cv_kj.order + 1, n_kj, i);
            CombinatorialVertex scoop_cv2(k, j, 2 * scoop_idx + cv_kj.order + 2, n_kj, i);
            interior_region_vertices.push_back(scoop_cv1);
            interior_region_vertices.push_back(scoop_cv2);
        }
        CombinatorialVertex cv_jk(j, k, smeared_curve_nc.scc.get(j, k, i), smeared_curve_nc.edge_ints[{j, k}], i);
        interior_region_vertices.push_back(cv_jk);
        return {CombFace(interior_region_vertices, CombFaceType::INTERIOR_CORNER_TYPE)};
    }
}

inline TriangleSoup diagonal_surface_construction(const array<Vec3, 4>& tet_positions,
                                                  const array<vector<double>, 6>& global_edge_isect_ts,
                                                  const CombFace& diagonal_comb_curve, bool scoop_mid_vertices,
                                                  double scoop_mid_vertex_bulge) {
    CombFace spiral_comb_face = diagonal_comb_curve;
    spiral_comb_face.type = CombFaceType::SPIRAL;
    if (scoop_mid_vertices) spiral_comb_face = add_mid_scoop_vertices_spiral(spiral_comb_face);
    return detail::emit_comb_face(spiral_comb_face, tet_positions, global_edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge);
}

inline TriangleSoup nondiagonal_surface_construction(const array<Vec3, 4>& tet_positions,
                                                     const array<vector<double>, 6>& global_edge_isect_ts,
                                                     const CombFace& nondiagonal_comb_curve, bool scoop_mid_vertices,
                                                     double scoop_mid_vertex_bulge) {
    auto [local_edge_isect_ts, curve_edge_occupations] =
        detail::extract_curve_edge_isect_ts(nondiagonal_comb_curve.vertices, global_edge_isect_ts);
    EdgeInts smeared_curve_edge_ints = combface_to_edgeints(nondiagonal_comb_curve);
    int interior_vertex = detail::get_interior_vertex(separated_vertices(smeared_curve_edge_ints));
    assert(interior_vertex >= 0 ||
           (smeared_curve_edge_ints[{0, 1}] % 2 == 0 && smeared_curve_edge_ints[{0, 2}] % 2 == 0 &&
            smeared_curve_edge_ints[{0, 3}] % 2 == 0 && smeared_curve_edge_ints[{1, 2}] % 2 == 0 &&
            smeared_curve_edge_ints[{1, 3}] % 2 == 0 && smeared_curve_edge_ints[{2, 3}] % 2 == 0));

    NormalCoordinates smeared_curve_nc(smeared_curve_edge_ints);
    if (smeared_curve_nc.triangle_ineq || !smeared_curve_nc.even_sum)
        throw std::invalid_argument("smeared_surface_construction: should satisfy even sum and NOT satisfy triangle inequality");

    vector<CombFace> comb_faces;
    auto append = [&comb_faces](const vector<CombFace>& section) { comb_faces.insert(comb_faces.end(), section.begin(), section.end()); };
    append(smeared_quads_construction(smeared_curve_nc, interior_vertex));
    append(interior_region_construction(smeared_curve_nc, interior_vertex));
    if (scoop_mid_vertices) append(tunnel_construction(smeared_curve_nc, interior_vertex));

    TriangleSoup curve_soup;
    for (CombFace face : comb_faces) {
        if (scoop_mid_vertices) face = add_mid_scoop_vertices(face, smeared_curve_nc);
        if (face.vertices.size() < 3) continue;
        curve_soup.add_local_soup(detail::emit_comb_face(face, tet_positions, local_edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge));
    }
    if (TriangleSoup::COMB_MERGE) curve_soup.remap_orders(global_edge_isect_ts, curve_edge_occupations, true);
    return curve_soup;
}

inline TriangleSoup evensum_surface_construction(const array<Vec3, 4>& tet_positions,
                                                 const array<vector<double>, 6>& global_edge_isect_ts,
                                                 const vector<CombFace>& non_normal_comb_curves,
                                                 bool scoop_mid_vertices, double scoop_mid_vertex_bulge) {
    TriangleSoup soup;
    for (auto& non_normal_comb_curve : non_normal_comb_curves) {
        EdgeInts evensum_edge_ints = combface_to_edgeints(non_normal_comb_curve);
        vector<int> separated_vs = separated_vertices(evensum_edge_ints);
        if (separated_vs.size() == 2) {
            soup.add_local_soup(diagonal_surface_construction(tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge));
        } else {
            soup.add_local_soup(nondiagonal_surface_construction(tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge));
        }
    }
    return soup;
}

// ============================================================================
// normal reconstruction  (from subgrid_MT/normal_reconstruction.cpp)
// ============================================================================
// Construction policy (NOT a runtime option): separate the corner disks from the
// recursive subdivision. This is the shipped behavior and should only be changed
// by a developer who understands the almost-normal recursion. It is deliberately
// a file-local constant rather than a plumbed-through flag.
static constexpr bool no_corner_recursion = true;

// forward declaration (mutual recursion normal_surface_construction <-> recursive_construction)
inline TriangleSoup normal_surface_construction(const array<Vec3, 4>& tet_positions,
                                                const array<vector<double>, 6>& normal_edge_isect_ts,
                                                bool verbose);

namespace detail {

inline void shift_orders_past_corners(TriangleSoup& soup, const NormalCoordinates& NC) {
    for (auto& sig : soup.signatures) {
        if (sig.type == CombVertexSigType::STANDALONE) continue;
        if (sig.i == 4 || sig.j == 4) continue;
        sig.order += NC.cc[sig.i];
        sig.edge_int = NC.edge_ints[{sig.i, sig.j}];
    }
}

}  // namespace detail

inline TriangleSoup corner_construction(const array<Vec3, 4>& tet_positions,
                                        const array<vector<double>, 6>& normal_edge_isect_ts,
                                        const NormalCoordinates& NC) {
    TriangleSoup soup;
    for (int i = 0; i < 4; i++) {
        vector<int> others = complement({i});
        int j = others[0], k = others[1], l = others[2];
        for (int order = 1; order <= NC.cc[i]; order++) {
            vector<CombinatorialVertex> triangle_comb_verts{
                CombinatorialVertex(i, j, order, NC.edge_ints[{i, j}]),
                CombinatorialVertex(i, k, order, NC.edge_ints[{i, k}]),
                CombinatorialVertex(i, l, order, NC.edge_ints[{i, l}])};
            CombFace comb_face(triangle_comb_verts, CombFaceType::TRIANGLE);
            soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
        }
    }
    return soup;
}

inline TriangleSoup diagonal_construction(const array<Vec3, 4>& tet_positions,
                                          const array<vector<double>, 6>& normal_edge_isect_ts,
                                          const NormalCoordinates& NC) {
    TriangleSoup soup;
    int b1 = 0;
    int b2 = NC.dc[{b1, 1}] != 0 ? 1 : NC.dc[{b1, 2}] != 0 ? 2 : NC.dc[{b1, 3}] != 0 ? 3 : -1;
    if (b2 == -1) return soup;
    vector<int> head_vs = complement({b1, b2});
    int h1 = head_vs[0], h2 = head_vs[1];
    int diag_count = NC.dc[{b1, b2}];
    int b1_cc = NC.cc[b1], b2_cc = NC.cc[b2];
    for (int t = 0; t < diag_count; t++) {
        CombinatorialVertex cv0(b1, h1, t + b1_cc + 1, NC.edge_ints[{b1, h1}]);
        CombinatorialVertex cv1(b1, h2, t + b1_cc + 1, NC.edge_ints[{b1, h2}]);
        CombinatorialVertex cv2(b2, h2, t + b2_cc + 1, NC.edge_ints[{b2, h2}]);
        CombinatorialVertex cv3(b2, h1, t + b2_cc + 1, NC.edge_ints[{b2, h1}]);
        CombFace comb_face({cv0, cv1, cv2, cv3}, CombFaceType::QUAD02);
        soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    }
    return soup;
}

inline TriangleSoup octagon_construction(const array<Vec3, 4>& tet_positions,
                                         const array<vector<double>, 6>& normal_edge_isect_ts,
                                         const NormalCoordinates& NC) {
    TriangleSoup soup;
    int b1 = 0;
    vector<int> other_vs{1, 2, 3};
    int vd1 = -1, vd2 = -1, vdd = -1;
    int d1 = -1, d2 = -1;
    for (int other_b1 : other_vs) {
        int tmp_d = NC.dc[{b1, other_b1}];
        if (tmp_d != 0) {
            if (vd1 == -1) { vd1 = other_b1; d1 = tmp_d; continue; }
            else { vd2 = other_b1; d2 = tmp_d; }
        } else vdd = other_b1;
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
    for (int t = 0; t < d; t++) {
        DerivedVertex inner_vertex_t(vd1_vd2_mid_vertex, b1_vdd_mid_vertex, (double)(t + 1) / (double)(d + 1));
        CombinatorialVertex b1vdd_1(b1, vdd, t + 1 + c_b1, NC.edge_ints[{b1, vdd}]);
        CombinatorialVertex b1vdd_2(vdd, b1, t + 1 + c_vdd, NC.edge_ints[{b1, vdd}]);
        CombinatorialVertex d1d2_1(vd1, vd2, d - t + c_vd1, NC.edge_ints[{vd1, vd2}]);
        CombinatorialVertex d1d2_2(vd2, vd1, d - t + c_vd2, NC.edge_ints[{vd1, vd2}]);
        CombinatorialVertex b1vd1(b1, vd1, t + 1 + c_b1, NC.edge_ints[{b1, vd1}]);
        CombinatorialVertex b1vd2(b1, vd2, t + 1 + c_b1, NC.edge_ints[{b1, vd2}]);
        CombinatorialVertex d1vdd(vd1, vdd, d - t + c_vd1, NC.edge_ints[{vd1, vdd}]);
        CombinatorialVertex d2vdd(vd2, vdd, d - t + c_vd2, NC.edge_ints[{vd2, vdd}]);
        CombFace octagon_comb_face({b1vdd_1, b1vd1, d1d2_1, d1vdd, b1vdd_2, d2vdd, d1d2_2, b1vd2}, CombFaceType::OCTAGON);
        octagon_comb_face.derived_vertex = inner_vertex_t;
        soup.add_local_soup(embed_comb_face(octagon_comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    }
    return soup;
}

inline TriangleSoup spiral_construction(const array<Vec3, 4>& tet_positions,
                                        const array<vector<double>, 6>& normal_edge_isect_ts,
                                        const NormalCoordinates& NC) {
    assert(NC.even_sum && NC.triangle_ineq);
    int i = 0, j = 1;
    CombinatorialVertex start_cv(i, j, NC.cc[i] + 1, NC.edge_ints[{i, j}]);
    vector<CombinatorialVertex> spiral_comb_verts;
    array<size_t, 4> dummy_global_inds = {0, 1, 2, 3};
    EdgeOccupations dummy_edge_occ;
    trace_generic_curve(start_cv, NC, dummy_global_inds, dummy_edge_occ, spiral_comb_verts, false);
    CombFace spiral_comb_face(spiral_comb_verts, CombFaceType::SPIRAL);
    TriangleSoup soup;
    soup.add_local_soup(embed_comb_face(spiral_comb_face, tet_positions, normal_edge_isect_ts, false, 0.0));
    return soup;
}

inline pair<Vec3, array<double, 4>>
compute_recursion_center_and_isect_ts(const array<Vec3, 4>& tet_positions,
                                      const array<vector<double>, 6>& normal_edge_isect_ts,
                                      const NormalCoordinates& NC) {
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
    Vec3 p01 = interpolate_comb_vertex(cv01, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p10 = interpolate_comb_vertex(cv10, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p02 = interpolate_comb_vertex(cv02, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p20 = interpolate_comb_vertex(cv20, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p03 = interpolate_comb_vertex(cv03, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p30 = interpolate_comb_vertex(cv30, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p12 = interpolate_comb_vertex(cv12, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p21 = interpolate_comb_vertex(cv21, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p13 = interpolate_comb_vertex(cv13, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p31 = interpolate_comb_vertex(cv31, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p23 = interpolate_comb_vertex(cv23, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 p32 = interpolate_comb_vertex(cv32, tet_positions, normal_edge_isect_ts, -1.0);
    Vec3 subd_center = (p01 + p10 + p02 + p20 + p03 + p30 + p12 + p21 + p13 + p31 + p23 + p32) / 12.0;
    double ts0_limit, ts1_limit, ts2_limit, ts3_limit;
    bool hit0 = segment_triangle_intersection(subd_center, tet_positions[0], p01, p02, p03, ts0_limit);
    bool hit1 = segment_triangle_intersection(subd_center, tet_positions[1], p10, p12, p13, ts1_limit);
    bool hit2 = segment_triangle_intersection(subd_center, tet_positions[2], p20, p21, p23, ts2_limit);
    bool hit3 = segment_triangle_intersection(subd_center, tet_positions[3], p30, p31, p32, ts3_limit);
    if (!hit0 || !hit1 || !hit2 || !hit3) {
        std::cout << "Warning: compute_recursion_center_and_isect_ts: segment-triangle intersection failed; falling back to greedy recursion intersections" << std::endl;
        ts0_limit = ts1_limit = ts2_limit = ts3_limit = 1.0;
    }
    return std::make_pair(subd_center, array<double, 4>{{ts0_limit, ts1_limit, ts2_limit, ts3_limit}});
}

inline TriangleSoup recursive_construction(const array<Vec3, 4>& tet_positions,
                                           const array<vector<double>, 6>& normal_edge_isect_ts,
                                           const NormalCoordinates& NC, bool verbose) {
    auto [d1, d2] = NC.get_d1_d2();
    if (d1 < d2) std::swap(d1, d2);
    array<int, 4> interior_edge_ints = {0, 0, 0, 0};
    int b1 = 0;
    interior_edge_ints[b1] = d1 - d2 + (no_corner_recursion ? 0 : NC.cc[b1]);
    array<int, 3> others_vs = {1, 2, 3};
    for (int other_v : others_vs) {
        if (NC.dc[{b1, other_v}] == 0) interior_edge_ints[other_v] = 2 * d2 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else if (NC.dc[{b1, other_v}] == d1) interior_edge_ints[other_v] = d1 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else if (NC.dc[{b1, other_v}] == d2) interior_edge_ints[other_v] = d2 + (no_corner_recursion ? 0 : NC.cc[other_v]);
        else throw std::logic_error("recursive_construction: diagonal cuts cant take any other values. NC values are off\n");
    }
    Vec3 new_interior_tet_vertex = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0;
    array<double, 4> interior_subd_ts_limit = {1.0, 1.0, 1.0, 1.0};
    if (no_corner_recursion) {
        auto [recursion_center, interior_isect_limits] = compute_recursion_center_and_isect_ts(tet_positions, normal_edge_isect_ts, NC);
        new_interior_tet_vertex = recursion_center;
        for (int i = 0; i < 4; i++) interior_subd_ts_limit[i] = interior_isect_limits[i];
    }
    TriangleSoup soup;
    for (array<int, 3> triplet : ALL_TET_TRIPLETS) {
        int i = triplet[0], j = triplet[1], k = triplet[2];
        array<Vec3, 4> interior_tet_positions = {new_interior_tet_vertex, tet_positions[i], tet_positions[j], tet_positions[k]};
        array<vector<double>, 6> interior_tet_edge_isect_ts;
        interior_tet_edge_isect_ts[0] = generate_uniform_isect_ts_single_edge(interior_edge_ints[i], interior_subd_ts_limit[i]);
        interior_tet_edge_isect_ts[1] = generate_uniform_isect_ts_single_edge(interior_edge_ints[j], interior_subd_ts_limit[j]);
        interior_tet_edge_isect_ts[2] = generate_uniform_isect_ts_single_edge(interior_edge_ints[k], interior_subd_ts_limit[k]);
        if (no_corner_recursion) {
            auto& ts_ij = normal_edge_isect_ts[edge_pair_to_index(i, j)];
            auto& ts_ik = normal_edge_isect_ts[edge_pair_to_index(i, k)];
            auto& ts_jk = normal_edge_isect_ts[edge_pair_to_index(j, k)];
            interior_tet_edge_isect_ts[3] = vector<double>(ts_ij.begin() + NC.cc[i], ts_ij.end() - NC.cc[j]);
            interior_tet_edge_isect_ts[4] = vector<double>(ts_ik.begin() + NC.cc[i], ts_ik.end() - NC.cc[k]);
            interior_tet_edge_isect_ts[5] = vector<double>(ts_jk.begin() + NC.cc[j], ts_jk.end() - NC.cc[k]);
        } else {
            interior_tet_edge_isect_ts[3] = normal_edge_isect_ts[edge_pair_to_index(i, j)];
            interior_tet_edge_isect_ts[4] = normal_edge_isect_ts[edge_pair_to_index(i, k)];
            interior_tet_edge_isect_ts[5] = normal_edge_isect_ts[edge_pair_to_index(j, k)];
        }
        auto sub_soup = normal_surface_construction(interior_tet_positions, interior_tet_edge_isect_ts, verbose);
        if (TriangleSoup::COMB_MERGE) {
            sub_soup.remap_vertex_indices({4, i, j, k});
            if (no_corner_recursion) detail::shift_orders_past_corners(sub_soup, NC);
        }
        soup.add_local_soup(sub_soup);
    }
    if (TriangleSoup::COMB_MERGE) soup.make_vertex_standalone(4);
    return soup;
}

inline TriangleSoup normal_surface_construction(const array<Vec3, 4>& tet_positions,
                                                const array<vector<double>, 6>& normal_edge_isect_ts,
                                                bool verbose) {
    EdgeInts normal_edge_ints(normal_edge_isect_ts);
    NormalCoordinates NC(normal_edge_ints);
    if (!NC.triangle_ineq || !NC.even_sum)
        throw std::runtime_error("Invalid normal coordinates: triangle inequality or even sum condition failed.");
    auto [d1, d2] = NC.get_d1_d2();
    TriangleSoup soup;
    if (d1 == 0 || d2 == 0) {
        auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
        auto diagonal_soup = diagonal_construction(tet_positions, normal_edge_isect_ts, NC);
        soup.add_local_soup(corner_soup);
        soup.add_local_soup(diagonal_soup);
    } else {
        if (d1 == d2 || std::gcd(d1, d2) == 1) {
            auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
            soup.add_local_soup(corner_soup);
            if (d1 == d2) {
                auto octagon_soup = octagon_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(octagon_soup);
            } else {
                auto spiral_soup = spiral_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(spiral_soup);
            }
        } else {
            if (no_corner_recursion) {
                auto corner_soup = corner_construction(tet_positions, normal_edge_isect_ts, NC);
                soup.add_local_soup(corner_soup);
            }
            auto recursive_soup = recursive_construction(tet_positions, normal_edge_isect_ts, NC, verbose);
            soup.add_local_soup(recursive_soup);
        }
    }
    return soup;
}

// ============================================================================
// primal entry  (from subgrid_MT/primal_reconstruction.cpp)
// ============================================================================
inline array<vector<double>, 6> unoccupied_edge_isect_ts(const array<vector<double>, 6>& edge_isect_ts,
                                                         const EdgeOccupations& edge_occupations) {
    array<vector<double>, 6> result;
    for (size_t eIdx = 0; eIdx < 6; ++eIdx) {
        int i = ALL_TET_PAIRS[eIdx].first, j = ALL_TET_PAIRS[eIdx].second;
        int n_ij = edge_isect_ts[eIdx].size();
        for (int k = 0; k < n_ij; ++k)
            if (!edge_occupations.get(i, j, k)) result[eIdx].push_back(edge_isect_ts[eIdx][k]);
    }
    return result;
}

inline TriangleSoup subgrid_surface(const array<Vec3, 4>& tet_positions, const array<size_t, 4>& tet_global_indices,
                                    const array<vector<double>, 6>& edge_isect_ts,
                                    const vector<CombFace>& open_boundary_curves,
                                    const vector<CombFace>& scoop_boundary_curves,
                                    const EdgeOccupations& non_normal_edge_occupations, bool scoop_mid_vertices,
                                    double scoop_mid_vertex_bulge, bool verbose = false) {
    TriangleSoup soup;
    // open curves are discarded
    TriangleSoup evensum_soup = evensum_surface_construction(tet_positions, edge_isect_ts, scoop_boundary_curves,
                                                             scoop_mid_vertices, scoop_mid_vertex_bulge);
    soup.add_local_soup(evensum_soup);
    auto normal_edge_isect_ts = unoccupied_edge_isect_ts(edge_isect_ts, non_normal_edge_occupations);
    auto normal_soup = normal_surface_construction(tet_positions, normal_edge_isect_ts, verbose);
    if (TriangleSoup::COMB_MERGE) normal_soup.remap_orders(edge_isect_ts, non_normal_edge_occupations, false);
    soup.add_local_soup(normal_soup);
    if (TriangleSoup::COMB_MERGE) {
        array<int, 4> idx_map = {(int)tet_global_indices[0], (int)tet_global_indices[1],
                                 (int)tet_global_indices[2], (int)tet_global_indices[3]};
        soup.remap_vertex_indices(idx_map);
        soup.make_type_standalone(CombVertexSigType::SCOOP_INTERIOR_STEINER);
    }
    return soup;
}

// ============================================================================
// Public output structs + folded entry points
// ============================================================================
struct PrimalResult {
    vector<Vec3> vertices;
    vector<vector<size_t>> faces;              // polygons (n-gons), indices into vertices
    vector<CombVertexKey> signatures;          // parallel to vertices; use signature_key() to merge
    bool non_even = false;                     // true if this tet is NOT even-sum (open curves present)
};

// Folds boundary_comb_curves inside; input matches construction_policy.md "Input".
inline PrimalResult subgrid_primal(const array<Vec3, 4>& tet_positions,
                                   const array<size_t, 4>& tet_global_indices,
                                   const array<vector<double>, 6>& edge_isect_ts,
                                   bool scoop_mid_vertices = true, double scoop_mid_vertex_bulge = 1e-4) {
    TriangleSoup::COMB_MERGE = true;
    array<int, 6> counts{};
    for (int e = 0; e < 6; e++) counts[e] = (int)edge_isect_ts[e].size();

    EdgeOccupations non_normal_edge_occupations;
    auto [open_curves, scoop_curves, normal_curves] =
        boundary_comb_curves(tet_global_indices, counts, non_normal_edge_occupations, /*populate_normal_curves=*/false);

    TriangleSoup soup = subgrid_surface(tet_positions, tet_global_indices, edge_isect_ts, open_curves, scoop_curves,
                                        non_normal_edge_occupations, scoop_mid_vertices, scoop_mid_vertex_bulge,
                                        /*verbose=*/false);

    PrimalResult out;
    out.vertices = std::move(soup.vertices);
    out.faces = std::move(soup.faces);
    out.signatures = std::move(soup.signatures);
    out.non_even = !open_curves.empty();
    return out;
}

// ============================================================================
// Dual: QEF solver  (from subgrid_MT/dual/qef_solver.cpp; Eigen ldlt -> hand-rolled 3x3)
// ============================================================================
namespace detail {

// Solve A x = b for a 3x3 matrix via the adjugate/determinant (analytic inverse).
// A is symmetric positive-definite here (sum of n n^T + reg*I, reg > 0), so this
// stands in for Eigen's A.ldlt().solve(b).
inline Vec3 solve3x3(const double A[3][3], const Vec3& b) {
    double c00 = A[1][1] * A[2][2] - A[1][2] * A[2][1];
    double c01 = A[1][0] * A[2][2] - A[1][2] * A[2][0];
    double c02 = A[1][0] * A[2][1] - A[1][1] * A[2][0];
    double det = A[0][0] * c00 - A[0][1] * c01 + A[0][2] * c02;
    // inverse (adjugate transposed) applied to b
    double i00 = c00 / det;
    double i01 = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det;
    double i02 = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
    double i10 = -c01 / det;
    double i11 = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
    double i12 = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det;
    double i20 = c02 / det;
    double i21 = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det;
    double i22 = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;
    return {i00 * b.x + i01 * b.y + i02 * b.z, i10 * b.x + i11 * b.y + i12 * b.z,
            i20 * b.x + i21 * b.y + i22 * b.z};
}

}  // namespace detail

inline vector<Vec3> QEF_from_boundary_polygons(const vector<vector<size_t>>& polygons,
                                               const vector<Vec3>& positions, const vector<Vec3>& normals,
                                               double regularizer_weight = 1e-3) {
    vector<Vec3> dual_points;
    for (const auto& poly : polygons) {
        if (poly.size() < 2) {
            dual_points.push_back(Vec3::undefined());
            continue;
        }
        Vec3 avg_pos = Vec3::zero();
        for (size_t vidx : poly) avg_pos += positions[vidx];
        avg_pos /= static_cast<double>(poly.size());
        double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        Vec3 b = Vec3::zero();
        for (size_t vidx : poly) {
            const Vec3& n = normals[vidx];
            const Vec3& p = positions[vidx];
            A[0][0] += n.x * n.x; A[0][1] += n.x * n.y; A[0][2] += n.x * n.z;
            A[1][0] += n.y * n.x; A[1][1] += n.y * n.y; A[1][2] += n.y * n.z;
            A[2][0] += n.z * n.x; A[2][1] += n.z * n.y; A[2][2] += n.z * n.z;
            double ndotp = dot(n, p);
            b += n * ndotp;
        }
        A[0][0] += regularizer_weight; A[1][1] += regularizer_weight; A[2][2] += regularizer_weight;
        b += regularizer_weight * avg_pos;
        dual_points.push_back(detail::solve3x3(A, b));
    }
    return dual_points;
}

inline vector<Vec3> get_centroids(const vector<vector<size_t>>& polygons, const vector<Vec3>& positions) {
    vector<Vec3> centroids;
    for (const auto& poly : polygons) {
        Vec3 centroid = Vec3::zero();
        for (size_t vidx : poly) centroid += positions[vidx];
        centroid /= static_cast<double>(poly.size());
        centroids.push_back(centroid);
    }
    return centroids;
}

inline Vec3 projected_dual_position(const Vec3& dual_p, const Vec3& centroid, const array<Vec3, 4>& tet_positions) {
    if (is_inside_tet(dual_p, tet_positions)) return dual_p;
    double eps = 1e-6;
    Vec3 tet_center = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0;
    Vec3 adjusted_centroid = centroid + 1e-6 * (tet_center - centroid);
    double isect_t;
    segment_tet_intersection(dual_p, adjusted_centroid, tet_positions, isect_t);
    Vec3 dir = adjusted_centroid - dual_p;
    return dual_p + (isect_t + eps) * dir;
}

inline vector<Vec3> projected_dual_positions(const vector<Vec3>& dual_positions, const vector<Vec3>& centroids,
                                             const array<Vec3, 4>& tet_positions) {
    vector<Vec3> projected_positions;
    for (size_t idx = 0; idx < dual_positions.size(); idx++) {
        Vec3 dp = dual_positions[idx];
        if (!dp.isDefined()) { projected_positions.push_back(dp); continue; }
        projected_positions.push_back(projected_dual_position(dp, centroids[idx], tet_positions));
    }
    return projected_positions;
}

inline vector<Vec3> dual_positions_from_isect_data(const array<Vec3, 4>& tet_positions,
                                                   const vector<vector<size_t>>& boundary_polygons,
                                                   const vector<Vec3>& isect_positions,
                                                   const vector<Vec3>& isect_normals, double regularizer_weight,
                                                   bool projected) {
    vector<Vec3> dual_positions = QEF_from_boundary_polygons(boundary_polygons, isect_positions, isect_normals, regularizer_weight);
    if (projected) {
        vector<Vec3> centroids = get_centroids(boundary_polygons, isect_positions);
        dual_positions = projected_dual_positions(dual_positions, centroids, tet_positions);
    }
    return dual_positions;
}

// ============================================================================
// Dual construction  (from subgrid_MT/dual/dual_construction.cpp)
// ============================================================================
inline vector<Vec3> normals_from_cv_polygons(const vector<CombFace>& curves,
                                             const array<vector<double>, 6>& edge_isect_ts,
                                             const array<vector<Vec3>, 6>& edge_isect_normals) {
    vector<Vec3> normals;
    for (const auto& curve : curves) {
        for (const auto& cv : curve.vertices) {
            int i = cv.i, j = cv.j;
            int edge_idx = edge_pair_to_index(i, j);
            size_t idx = (i < j) ? (cv.order - 1) : (cv.edge_int - cv.order);
            normals.push_back(edge_isect_normals[edge_idx][idx]);
        }
    }
    return normals;
}

inline TriangleSoup build_dual_local_soup(const vector<CombFace>& curves, const array<Vec3, 4>& tet_positions,
                                          const array<vector<double>, 6>& edge_isect_ts) {
    TriangleSoup soup;
    size_t vert_offset = 0;
    for (const auto& curve : curves) {
        size_t n = curve.vertices.size();
        vector<size_t> polygon(n);
        auto fpe = curve_segment_faces(curve);
        for (size_t vi = 0; vi < n; vi++) {
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

inline TriangleSoup dual_subgrid_surface(const array<Vec3, 4>& tet_positions, const array<size_t, 4>& tet_indices,
                                         const array<vector<double>, 6>& edge_isect_ts,
                                         const array<vector<Vec3>, 6>& edge_isect_normals,
                                         const vector<CombFace>& open_curves, const vector<CombFace>& scoop_curves,
                                         const vector<CombFace>& normal_curves, double reg_alpha, bool project_duals) {
    (void)open_curves;
    vector<CombFace> curves;
    for (const auto& c : scoop_curves) curves.push_back(c);
    for (const auto& c : normal_curves) curves.push_back(c);
    if (curves.empty()) return {};

    auto normals = normals_from_cv_polygons(curves, edge_isect_ts, edge_isect_normals);
    TriangleSoup soup = build_dual_local_soup(curves, tet_positions, edge_isect_ts);
    soup.dual_positions = dual_positions_from_isect_data(tet_positions, soup.faces, soup.vertices, normals, reg_alpha, project_duals);
    array<int, 4> idx_map = {(int)tet_indices[0], (int)tet_indices[1], (int)tet_indices[2], (int)tet_indices[3]};
    soup.remap_vertex_indices(idx_map);
    return soup;
}

// ============================================================================
// Dual entry (folded)
// ============================================================================
struct DualResult {
    vector<Vec3> vertices;
    vector<vector<size_t>> faces;                    // boundary loops (one polygon per dual point)
    vector<vector<array<size_t, 3>>> faces_per_edge; // parallel to faces; tet face {i,j,k} per polygon edge (global)
    vector<Vec3> dual_positions;                     // parallel to faces; one QEF dual point per polygon
    vector<CombVertexKey> signatures;                // parallel to vertices; all NORMAL type
    bool non_even = false;
    bool non_normal = false;
};

inline DualResult subgrid_dual(const array<Vec3, 4>& tet_positions, const array<size_t, 4>& tet_global_indices,
                               const array<vector<double>, 6>& edge_isect_ts,
                               const array<vector<Vec3>, 6>& edge_isect_normals, double reg_alpha = 0.1,
                               bool project_duals = false) {
    TriangleSoup::COMB_MERGE = true;
    array<int, 6> counts{};
    for (int e = 0; e < 6; e++) counts[e] = (int)edge_isect_ts[e].size();

    EdgeOccupations edge_occupations;
    auto [open_curves, scoop_curves, normal_curves] =
        boundary_comb_curves(tet_global_indices, counts, edge_occupations, /*populate_normal_curves=*/true);

    TriangleSoup soup = dual_subgrid_surface(tet_positions, tet_global_indices, edge_isect_ts, edge_isect_normals,
                                             open_curves, scoop_curves, normal_curves, reg_alpha, project_duals);
    DualResult out;
    out.vertices = std::move(soup.vertices);
    out.faces = std::move(soup.faces);
    out.faces_per_edge = std::move(soup.faces_per_edge);
    out.dual_positions = std::move(soup.dual_positions);
    out.signatures = std::move(soup.signatures);
    out.non_even = !open_curves.empty();
    out.non_normal = !(open_curves.empty() && scoop_curves.empty());
    return out;
}

}  // namespace subgrid_mt
