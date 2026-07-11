#pragma once

#include <array>
#include <map>
#include <utility>

#include "common/utils.h"

struct EdgeInts{
    // Access as ei[{i,j}] with unordered pair {i,j} where i!=j and i,j in [0,3].
    // e.g. ei[{0,1}] is e01
    std::map<std::pair<int,int>, int> data;

    static std::pair<int,int> key(int i, int j){ if (i>j) std::swap(i,j); return {i,j}; }

    
    // Constructors
    EdgeInts() = default;
    EdgeInts(int n_01, int n_02, int n_03, int n_12, int n_13, int n_23){
        set(0,1,n_01); set(0,2,n_02); set(0,3,n_03); set(1,2,n_12); set(1,3,n_13); set(2,3,n_23);
    }
    EdgeInts(std::array<int,6> counts){
        set(0,1,counts[0]); set(0,2,counts[1]); set(0,3,counts[2]);set(1,2,counts[3]); set(1,3,counts[4]); set(2,3,counts[5]);
    }
    EdgeInts(const std::array<std::vector<double>,6> &isects){
        set(0,1,isects[0].size()); set(0,2,isects[1].size()); set(0,3,isects[2].size());set(1,2,isects[3].size()); set(1,3,isects[4].size()); set(2,3,isects[5].size());
    }

    // convert to array
    std::array<int,6> to_array() const {
        std::array<int,6> counts;
        counts[0] = get(0,1); counts[1] = get(0,2); counts[2] = get(0,3);
        counts[3] = get(1,2); counts[4] = get(1,3); counts[5] = get(2,3);
        return counts;
    }

    // primitives
    int get(int i, int j) const {
        auto it = data.find(key(i,j));
        return it == data.end() ? 0 : it->second;
    }

    void set(int i, int j, size_t v){ data[key(i,j)] = v; }

    // overload []
    int operator[](std::pair<int,int> p) const {
        return get(p.first, p.second);
    }
    int& operator[](std::pair<int,int> p) {
        auto k = key(p.first, p.second);
        return data[k]; // inserts 0 if absent
    }

    // In-place subtraction
    EdgeInts& operator-=(const EdgeInts& other) {
        for (const auto& kv : other.data) {
            auto& v = data[kv.first]; // inserts 0 if absent
            v -= kv.second;
        }
        return *this;
    }
    // Value-returning subtraction
    friend EdgeInts operator-(EdgeInts lhs, const EdgeInts& rhs) {
        lhs -= rhs; // reuse operator-=
        return lhs;
    }

    // In-place Addition
    EdgeInts& operator+=(const EdgeInts& other) {
        for (const auto& kv : other.data) {
            auto& v = data[kv.first]; // inserts 0 if absent
            v += kv.second;
        }
        return *this;
    }
    // Value-returning addition
    friend EdgeInts operator+(EdgeInts lhs, const EdgeInts& rhs) {
        lhs += rhs; // reuse operator+=
        return lhs;
    }

    // log
    void print() const {
        for (const auto& kv : data) {
            std::cout << "e{" << kv.first.first << "," << kv.first.second << "}: " << kv.second << "  ,  ";
        }
        std::cout << std::endl;
    }
};


struct SemiCornerCuts{
    // Access as sc[i][{j,k}] with i in [0,3], unordered pair {j,k} where j!=k and j,k in [0,3] and j!=i,k!=i.
    std::array<std::map<std::pair<int,int>, int>, 4> data;

    static std::pair<int,int> key(int j, int k){ if (j>k) std::swap(j,k); return {j,k}; }

    int get(int i, int j, int k) const {
        auto it = data[i].find(key(j,k));
        return it == data[i].end() ? 0 : it->second;
    }

    void set(int i, int j, int k, int v){ data[i][key(j,k)] = v; }

    // overload []
    const std::map<std::pair<int,int>, int>& operator[](size_t i) const { return data[i]; }
    std::map<std::pair<int,int>, int>& operator[](size_t i) { return data[i]; }

    // log; only non-zero entries
    void print() const {
        for (int i = 0; i < 4; i++){
            for (const auto& kv : data[i]) {
                if (kv.second != 0){
                    std::cout << "scc{" << i << "," << kv.first.first << "," << kv.first.second << "}: " << kv.second << "  ,  ";
                }
            }
        }
        std::cout << std::endl;
    }
};

struct CornerCuts{
    int data[4];

    int get(int i) const { return data[i]; }
    void set(int i, int v) { data[i] = v; }

    // overload []
    int operator[](int i) const { return data[i]; }
    int& operator[](int i) { return data[i]; }

    // log; only non-zero entries
    void print() const {
        for (int i = 0; i < 4; i++){
            if (data[i] != 0){
                std::cout << "cc{" << i << "}: " << data[i] << "  ,  ";
            }
        }
        std::cout << std::endl;
    }
};

struct DiagonalCuts{
    // access as dc[{i,j}] where the order does not matter, also dc[{i,j}] = dc[{k,l}] where {i,j} and {k,l} are complementary pairs
    std::map<std::pair<int,int>, int> data;

    static std::pair<int,int> key(int i, int j){ if (i>j) std::swap(i,j); return {i,j}; }

    int get(int i, int j) const {
        auto it = data.find(key(i,j));
        return it == data.end() ? 0 : it->second;
    }

    void set(int i, int j, int v){ 
        data[key(i,j)] = v;
        // set complementary pair too
        for (int k = 0; k < 4; k++){
            if (k == i || k == j) continue;
            for (int l = k+1; l < 4; l++){
                if (l == i || l == j) continue;
                data[key(k,l)] = v;
            }
        }
    }

    // overload []
    int operator[](std::pair<int,int> p) const { return get(p.first, p.second); }
    int& operator[](std::pair<int,int> p) { return data[key(p.first, p.second)]; }

    // log; only non-zero entries
    void print() const {
        for (const auto& kv : data) {
            if (kv.second != 0){
                std::cout << "dc{" << kv.first.first << "," << kv.first.second << "}: " << kv.second << "  ,  ";
            }
        }
        std::cout << std::endl;
    }
};

// represents a combinatorial vertex in a single tet with indices 0, 1, 2, 3
// $$$ same as InterpolatedPoint but for topological construction phase $$$ 
struct CombinatorialVertex{ 
    int i = -1, j = -1, order = -1, edge_int = 0; // vertex is on edge {i,j} at position 'order' along the edge
    int k = -1; // the other vertex in the face containing this vertex
    int l = -1; // also use k as a tet index

    bool scoop_steiner = false; // when true, then k is the face side
    bool scoop_interior_steiner = false; // when true both k and l are used for shifting
    
    // Constructors
    CombinatorialVertex() = default;
    CombinatorialVertex(int i, int j, int order, int edge_int)
    : i(i), j(j), order(order), edge_int(edge_int){}
    CombinatorialVertex(int i, int j, int order, int edge_int, int k) // only used for tracing scoops
    : i(i), j(j), order(order), edge_int(edge_int), k(k) {}
    
    CombinatorialVertex(int i, int j, int order, int edge_int, int k, int l) // only used building the twin map
    : i(i), j(j), order(order), edge_int(edge_int), k(k), l(l) {}
    
    // // for handling auxilary inner vertices!!!
    // // TODO: handle this better!!!!!
    // std::array<double, 4> ts = {0.0, 0.0, 0.0, 0.0};
    // // Constructors for the bad inner vertices!!
    // CombinatorialVertex(int i, int j, int k, int l, std::array<double, 4> ts_)
    // : i(i), j(j), k(k), l(l), ts(ts_) {}

    bool operator==(const CombinatorialVertex& other) const {
        if ( (i == other.i) && (j == other.j) && (order == other.order))
            return true;
        if ( (i == other.j) && (j == other.i) && (order == edge_int - other.order + 1))
            return true;
        return false;
    }

    bool same_edge_as(const CombinatorialVertex& other) const {
        return ( (i == other.i && j == other.j) || (i == other.j && j == other.i) );
    }

    // log
    void print() const {
        std::cout << "CombinatorialVertex(i=" << i << ", j=" << j << ", order=" << order << ", edge_int=" << edge_int << ", k=" << k << ")";
    }
};


// derived vertex:
// represents a vertex that is derived from some first order combinatorial vertices
// Examples:
// -- A vertex created by averaging two combinatorial vertices on the same edge after a scoop
// -- Steiner vertex created by averaging a single spiral of combinatorial vertices
// -- Ordered Steiner vertices of normal octagons
struct DerivedVertex{
    std::vector<CombinatorialVertex> parents; // the combinatorial vertices that define this derived vertex
    std::vector<double> weights; // the weights for averaging the parents, should sum to 1

    // Constructors
    DerivedVertex() = default;
    // equal weights by default
    DerivedVertex(const std::vector<CombinatorialVertex>& parents) 
    : parents(parents), weights(std::vector<double>((double)parents.size(), 1.0 / (double)parents.size())) {}
    // custom weights
    DerivedVertex(const std::vector<CombinatorialVertex>& parents, const std::vector<double>& weights)
    : parents(parents), weights(weights) {}
    
    // from two other derived vertices by weights v1 * (1 - t1) + v2 * t1; merge the parents and weights
    DerivedVertex(const DerivedVertex& dv1, const DerivedVertex& dv2, double t1) {
        parents = dv1.parents;
        parents.insert(parents.end(), dv2.parents.begin(), dv2.parents.end());
        
        auto weights1 = dv1.weights;
        for (double& w1: weights1){
            w1 *= (1. - t1);
        }
        weights = weights1;

        auto weights2 = dv2.weights;
        for (double& w2: weights2){
            w2 *= t1;
        }
        weights.insert(weights.end(), weights2.begin(), weights2.end());
    }
};

// combindatorial face is a vector of combinatorial vertices 
// an enum for their type; this determines how they are converted to actual geometry later
enum class CombFaceType {
    OPEN, CLOSED_NON_NORMAL, CLOSED_NORMAL,
    TRIANGLE, QUAD02, QUAD13, 
    OCTAGON, SPIRAL, 
    // non-normal types
    SMEARED_QUAD, INTERIOR_HEXAGON, INTERIOR_CORNER_TYPE, 
    // CORNER_TYPE_TRIANGLE, // merged into INTERIOR_CORNER_TYPE
    TUNNEL_QUAD
};
struct CombFace{
    std::vector<CombinatorialVertex> vertices;
    CombFaceType type;
    
    // derived vertex for steiner point addition
    DerivedVertex derived_vertex; 

    // Constructors
    CombFace() = default;
    CombFace(const std::vector<CombinatorialVertex>& vertices, CombFaceType type)
    : vertices(vertices), type(type) {}

    // log
    void print() const {
        std::cout << "CombFace(type=";
        switch(type){
            case CombFaceType::TRIANGLE: std::cout << "TRIANGLE"; break;
            case CombFaceType::QUAD02: std::cout << "QUAD02"; break;
            case CombFaceType::QUAD13: std::cout << "QUAD13"; break;
            case CombFaceType::SPIRAL: std::cout << "SPIRAL"; break;
            case CombFaceType::OCTAGON: std::cout << "OCTAGON"; break;
            case CombFaceType::SMEARED_QUAD: std::cout << "SMEARED_QUAD"; break;
            case CombFaceType::INTERIOR_HEXAGON: std::cout << "INTERIOR_HEXAGON"; break;
            case CombFaceType::INTERIOR_CORNER_TYPE: std::cout << "INTERIOR_CORNER_TYPE"; break;
            case CombFaceType::TUNNEL_QUAD: std::cout << "TUNNEL_QUAD"; break;
            case CombFaceType::CLOSED_NON_NORMAL: std::cout << "CLOSED_NON_NORMAL"; break;
            case CombFaceType::CLOSED_NORMAL: std::cout << "CLOSED_NORMAL"; break;
            case CombFaceType::OPEN: std::cout << "OPEN"; break;
        }
        std::cout << ", vertices=[";
        for (const auto& v : vertices){
            v.print();
            std::cout << ", ";
        }
        std::cout << "])" << std::endl;
    }
};



struct EdgeOccupations{
    // Access as e[{i,j}] with unordered pair {i,j} where i!=j and i,j in [0,3].
    // each entry is a vector of booleans indicating occupation states along the edge
    // size of each entry is edge_ints[{i,j}] that should be given at some point
    std::map<std::pair<int,int>, std::vector<bool>> data;

    static std::pair<int,int> key(int i, int j){ if (i>j) std::swap(i,j); return {i,j}; }
    
    EdgeOccupations() = default;
    EdgeOccupations(const EdgeInts& edge_ints){
        // initialize data vectors to false
        for (const auto& kv : edge_ints.data){
            data[kv.first] = std::vector<bool>(kv.second, false);
        }
    }

    bool
    get(int i, int j, int index) const {
        if (i > j) { // i,j will be swapped, so index needs to be adjusted
            index = data.at(key(i,j)).size() - 1 - index;
        }
        return data.at(key(i,j))[index];
    }

    bool
    get(CombinatorialVertex cv) const {
        int i = cv.i;
        int j = cv.j;
        int index = cv.order - 1; // convert to 0-based
        if (i > j) { // i,j will be swapped, so index needs to be adjusted
            index = data.at(key(i,j)).size() - 1 - index;
        }
        return data.at(key(i,j))[index];
    }

    void 
    set(int i, int j, int index, bool value){
        if (i > j) { // i,j will be swapped, so index needs to be adjusted
            index = data[key(i,j)].size() - 1 - index;
        }
        data[key(i,j)][index] = value;
    }

    void 
    set(CombinatorialVertex cv, bool value){
        int i = cv.i;
        int j = cv.j;
        int index = cv.order - 1; // convert to 0-based
        if (i > j) { // i,j will be swapped, so index needs to be adjusted
            index = data[key(i,j)].size() - 1 - index;
        }
        data[key(i,j)][index] = value;
    }

    // get the number of occupied positions along edge {i,j} as EdgeInts
    EdgeInts to_edge_ints(bool count_occupied) const {
        EdgeInts result;
        for (const auto& kv : data){
            int count = 0;
            for (bool occupied : kv.second){
                if (occupied && count_occupied) count++;
                if (!occupied && !count_occupied) count++;
            }
            result[kv.first] = count;
        }
        return result;
    }



    // log
    void print() const {
        for (const auto& kv : data) {
            std::cout << "Edge {" << kv.first.first << "," << kv.first.second << "}: ";
            for (size_t idx = 0; idx < kv.second.size(); idx++) {
                std::cout << (kv.second[idx] ? "1" : "0");
                if (idx + 1 < kv.second.size()) std::cout << ",";
            }
            std::cout << "  ;  ";
        }
        std::cout << std::endl << std::endl;
    }
};


struct CVPool {
    // look into the pool and only create new if does not exist
    std::vector<CombinatorialVertex> cv_pool; 
    std::vector<std::tuple<int,int,int>> cv_keys;
    std::map<std::tuple<int,int,int>, size_t> cv_map; // i, j, order

    std::tuple<int, int, int> cv_key(const CombinatorialVertex& cv) {
        int i = cv.i, j = cv.j;
        int order = cv.order;
        if (i > j) {
            std::swap(i, j);  // normalize (i,j) ordering
            order = cv.edge_int - order + 1; // reverse order accordingly
        }
        return std::make_tuple(i, j, order);  // k is -1 if not a scoop vertex
    }

    // finder in pool
    size_t index(const CombinatorialVertex& cv) {
        auto key = cv_key(cv);
        // check if exists; otherwise add
        if (cv_map.find(key) == cv_map.end()){
            size_t new_idx = cv_pool.size();
            cv_pool.push_back(cv);
            cv_keys.push_back(key);
            cv_map[key] = new_idx;
            return new_idx;
        }
        else {
            return cv_map[key];
        }
    };

    // size of pool
    size_t size() { return cv_pool.size(); }

    // clear pool
    void clear() {
        cv_pool.clear();
        cv_keys.clear();
        cv_map.clear();
    }
};
