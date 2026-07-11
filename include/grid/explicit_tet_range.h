#pragma once

#include "common/utils.h"
#include "common/npz_reader.h"
#include <vector>
#include <stdexcept>

// Iterates over an explicit tet mesh loaded from an NPZ file.
// Same duck-typed interface as TetGridRange: begin()/end() yield objects
// with .indices (array<size_t,4>) and .positions (array<Vector3,4>).

class ExplicitTetRange {
public:
    struct TetData {
        std::array<size_t, 4> indices;
        std::array<Vector3, 4> positions;
    };

private:
    std::vector<TetData> tets;

public:
    ExplicitTetRange() = default;

    // Load from NPZ arrays: vertices (V,3) float64, tets (T,4) int32
    ExplicitTetRange(const NpyArray& vertices_arr, const NpyArray& tets_arr) {
        if (vertices_arr.shape.size() != 2 || vertices_arr.shape[1] != 3)
            throw std::runtime_error("ExplicitTetRange: vertices must be (V, 3)");
        if (tets_arr.shape.size() != 2 || tets_arr.shape[1] != 4)
            throw std::runtime_error("ExplicitTetRange: tets must be (T, 4)");

        size_t num_verts = vertices_arr.shape[0];
        size_t num_tets = tets_arr.shape[0];
        const double* verts = vertices_arr.as_float64();
        const int32_t* tet_indices = tets_arr.as_int32();

        tets.resize(num_tets);
        for (size_t t = 0; t < num_tets; ++t) {
            for (int i = 0; i < 4; ++i) {
                int32_t idx = tet_indices[t * 4 + i];
                if (idx < 0 || (size_t)idx >= num_verts)
                    throw std::runtime_error("ExplicitTetRange: tet index out of range");
                tets[t].indices[i] = (size_t)idx;
                tets[t].positions[i] = Vector3{
                    verts[idx * 3 + 0],
                    verts[idx * 3 + 1],
                    verts[idx * 3 + 2]
                };
            }
        }
    }

    using iterator = typename std::vector<TetData>::const_iterator;
    iterator begin() const { return tets.begin(); }
    iterator end() const { return tets.end(); }
    size_t total_tet_count() const { return tets.size(); }
};
