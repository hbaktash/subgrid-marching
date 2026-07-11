#pragma once

#include "nc/normalCoord_containers.h"

// corner and diagonal cuts to residual edge ints
EdgeInts corner_and_diagonal_cuts_to_edge_ints(
    const CornerCuts& cc,
    const DiagonalCuts& dc
);


class NormalCoordinates{
    public:
        EdgeInts edge_ints;
        
        // semicorner cuts
        SemiCornerCuts scc;
        // corner cuts
        CornerCuts cc;
        // diagonal cuts; not putting permutations in anymore
        DiagonalCuts dc;
        // scoops per edge, per triangle; scoop[i,j,k] is scoops on edge {j,k} in triangle opposite vertex i
        // Note: only one of scoop[i,j,k], scoop[j,i,k], scoop[k,i,j] is nonzero
        SemiCornerCuts scoop;

        bool triangle_ineq = true,
             even_sum = true;
        
        NormalCoordinates() = default;
        NormalCoordinates(EdgeInts &edge_ints)
        : edge_ints(edge_ints) {
            compute_NC();
        }
        NormalCoordinates(const EdgeInts &edge_ints)
        : edge_ints(edge_ints) {
            compute_NC();
        }
        
        // convert edge ints to NC data
        void compute_NC();
        // subroutines for compute_NC
        void compute_semi_corner_cuts_and_scoops();
        void compute_corner_cuts_and_reduced_semi_corner_cuts(SemiCornerCuts& scc_reduced);
        void compute_diagonal_cuts(const SemiCornerCuts& scc_reduced);

        // Scoop stuff
        EdgeInts scoop_reduction();

        // helper to trace scoop curves
        CombinatorialVertex find_next_vertex(CombinatorialVertex v0);
        // trace forever to fill occupancies
        void trace_forever(CombinatorialVertex v0, EdgeOccupations &edge_occupations);

        // get d1, d2
        std::pair<int, int> get_d1_d2() const;
};

