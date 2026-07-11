#pragma once

#include "common/utils.h"

// Lazy tet grid iterator: n³ cubes × 5 tets each, computed on-the-fly.
// dx = 1/(n-1) so vertices 0..n-1 span [-1,+1] exactly; vertex n extends
// slightly past +1. This avoids grid-surface coincidence for axis-aligned SDFs.
class TetGridIterator {
private:
    size_t n;  // grid resolution; creates n x n x n cubes
    double dx;
    bool center_and_normalize;
    
    // Current iteration state
    size_t ix, iy, iz;
    size_t tet_index;  // which tet within current cube (0-4)
    bool is_end;
    
    // Helper to compute a vertex position
    Vector3 compute_position(int x, int y, int z) const {
        Vector3 p{x * dx, y * dx, z * dx};
        if (center_and_normalize) {
            Vector3 center = Vector3::constant(0.5);
            p = (p - center) * 2.0;  // fit in [-1,1]^3
        }
        return p;
    }
    
    // Get the 8 corner indices of current cube
    // with the following layout:
    //              B----C
    //             /|   /|
    //            A----D |
    //            | F--|-G
    //            |/   |/
    //            E----H
    std::array<size_t, 8> get_cube_corner_indices() const {
        return {{
            iz * (n+1) * (n+1) + iy * (n+1) + ix,              // A
            iz * (n+1) * (n+1) + iy * (n+1) + (ix+1),          // B
            iz * (n+1) * (n+1) + (iy+1) * (n+1) + (ix+1),      // C
            iz * (n+1) * (n+1) + (iy+1) * (n+1) + ix,          // D
            (iz+1) * (n+1) * (n+1) + iy * (n+1) + ix,          // E
            (iz+1) * (n+1) * (n+1) + iy * (n+1) + (ix+1),      // F
            (iz+1) * (n+1) * (n+1) + (iy+1) * (n+1) + ix,      // G
            (iz+1) * (n+1) * (n+1) + (iy+1) * (n+1) + (ix+1)   // H
        }};
    }
    
    // Get the 5 tets for the current cube
    std::array<std::array<size_t, 4>, 5> get_cube_tets() const {
        auto [A, B, C, D, E, F, G, H] = get_cube_corner_indices();
        
        if ((ix + iy + iz) % 2 == 0) {
            return {{
                {A, C, F, G},  // ACFG
                {D, A, G, C},  // DACG
                {B, A, C, F},  // BACF
                {E, F, G, A},  // EFGA
                {H, F, C, G}   // HFCG
            }};
        } else {
            return {{
                {B, E, H, D},  // BEHD
                {B, E, A, D},  // BEAD
                {B, C, H, D},  // BCHD
                {B, E, H, F},  // BEHF
                {G, E, H, D}   // GEHD
            }};
        }
    }
    
    // Map corner index to world coordinates
    std::array<int, 3> index_to_coords(size_t idx) const {
        int z = idx / ((n+1) * (n+1));
        idx %= ((n+1) * (n+1));
        int y = idx / (n+1);
        int x = idx % (n+1);
        return {x, y, z};
    }
    
public:
    // Data structure returned by dereference
    struct TetData {
        std::array<size_t, 4> indices;
        std::array<Vector3, 4> positions;
    };
    
    TetGridIterator(size_t n_in, bool center_normalize, bool end_iter = false)
        : n(n_in), 
          dx(n_in > 1 ? 1.0 / static_cast<double>(n_in - 1) : 1.0),
          center_and_normalize(center_normalize),
          ix(0), iy(0), iz(0), tet_index(0), is_end(end_iter) 
    {}
    
    // Dereference operator: compute current tet data on-the-fly
    TetData operator*() const {
        auto tet_indices = get_cube_tets()[tet_index];
        std::array<Vector3, 4> tet_positions;
        
        for (int i = 0; i < 4; ++i) {
            auto [x, y, z] = index_to_coords(tet_indices[i]);
            tet_positions[i] = compute_position(x, y, z);
        }
        
        return {tet_indices, tet_positions};
    }
    
    // Prefix increment
    TetGridIterator& operator++() {
        if (is_end) return *this;
        
        tet_index++;
        if (tet_index >= 5) {
            tet_index = 0;
            ix++;
            if (ix >= n) {
                ix = 0;
                iy++;
                if (iy >= n) {
                    iy = 0;
                    iz++;
                    if (iz >= n) {
                        is_end = true;
                    }
                }
            }
        }
        return *this;
    }
    
    bool operator==(const TetGridIterator& other) const {
        return is_end == other.is_end &&
               (is_end || (ix == other.ix && iy == other.iy && iz == other.iz && tet_index == other.tet_index));
    }

    bool operator!=(const TetGridIterator& other) const {
        return !(*this == other);
    }

    size_t total_tet_count() const { return static_cast<size_t>(n) * n * n * 5; }
};

// Range wrapper for range-based for loops
class TetGridRange {
    size_t n;
    bool center_and_normalize;
    
public:
    TetGridRange(size_t n_in, bool cnorm = true) 
        : n(n_in), center_and_normalize(cnorm) {}
    
    TetGridIterator begin() const { 
        return TetGridIterator(n, center_and_normalize, false); 
    }
    
    TetGridIterator end() const { 
        return TetGridIterator(n, center_and_normalize, true); 
    }
    
    size_t total_tet_count() const {
        return n * n * n * 5;
    }
};