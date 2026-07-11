#pragma once

#include "common/utils.h"
#include "nc/normalCoord_containers.h"

#include <algorithm>
#include <array>
#include <random>
#include <set>
#include <vector>

// ── helpers ───────────────────────────────────────────────────────────────────

// Apply vertex permutation to an edge-coordinate 6-tuple.
// Under perm, the count on original edge {i,j} moves to new edge {perm[i],perm[j]}.
// Returns the 6-tuple in canonical edge order (01,02,03,12,13,23).
inline std::array<int,6> apply_perm(const std::array<int,6>& e,
                                     const std::array<int,4>& perm)
{
    std::array<int,6> result{};
    for (int idx = 0; idx < 6; ++idx) {
        int i = ALL_TET_PAIRS[idx].first;
        int j = ALL_TET_PAIRS[idx].second;
        int pi = perm[i], pj = perm[j];
        if (pi > pj) std::swap(pi, pj);
        result[edge_pair_to_index(pi, pj)] = e[idx];
    }
    return result;
}

// Canonical (lex-smallest) representative of an edge-coord tuple under all 24
// vertex permutations (S4 acts faithfully on the tet's 6 edges).
inline std::array<int,6> canonicalize(const std::array<int,6>& e) {
    std::array<int,4> perm = {0, 1, 2, 3};
    std::array<int,6> best = e;
    do {
        auto candidate = apply_perm(e, perm);
        if (candidate < best) best = candidate;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return best;
}

// ── generators ────────────────────────────────────────────────────────────────

// All distinct edge-coordinate 6-tuples with values in [0, max_k], up to
// tetrahedral (S4) symmetry.
inline std::vector<EdgeInts> enumerate_up_to_symmetry(int max_k) {
    std::set<std::array<int,6>> seen;
    const int base = max_k + 1;
    long long total = 1;
    for (int i = 0; i < 6; ++i) total *= base;

    for (long long n = 0; n < total; ++n) {
        std::array<int,6> e{};
        long long val = n;
        for (int i = 0; i < 6; ++i) { e[i] = (int)(val % base); val /= base; }
        seen.insert(canonicalize(e));
    }

    std::vector<EdgeInts> result;
    result.reserve(seen.size());
    for (const auto& arr : seen)
        result.emplace_back(arr);
    return result;
}

// n_samples random 6-tuples with values in [0, max_k], fixed seed.
// Discards samples where all edge counts are < min_k (since those are already
// covered by the exhaustive enumeration). Keeps drawing until n_samples accepted.
inline std::vector<EdgeInts> random_edge_coords(int max_k, int n_samples,
                                                 int min_k = 0,
                                                 unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, max_k);
    std::vector<EdgeInts> result;
    result.reserve(n_samples);
    while ((int)result.size() < n_samples) {
        std::array<int,6> e{};
        for (int i = 0; i < 6; ++i) e[i] = dist(rng);
        if (min_k > 0) {
            int mx = *std::max_element(e.begin(), e.end());
            if (mx < min_k) continue;
        }
        result.emplace_back(e);
    }
    return result;
}

// ── geometry generators ───────────────────────────────────────────────────────

// Generate sorted random t-values in (0,1) for each edge, count matching EdgeInts.
inline std::array<std::vector<double>,6>
random_isect_ts(const EdgeInts& e, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0.01, 0.99);
    std::array<std::vector<double>,6> ts;
    for (int idx = 0; idx < 6; ++idx) {
        int i = ALL_TET_PAIRS[idx].first;
        int j = ALL_TET_PAIRS[idx].second;
        int n = e.get(i, j);
        ts[idx].resize(n);
        for (int k = 0; k < n; ++k) ts[idx][k] = dist(rng);
        std::sort(ts[idx].begin(), ts[idx].end());
    }
    return ts;
}

// Fixed unit tet for single-tet tests.
inline std::array<Vector3,4> standard_tet_positions() {
    return {Vector3{0,0,0}, Vector3{1,0,0}, Vector3{0,1,0}, Vector3{0,0,1}};
}

// ── filter helpers ────────────────────────────────────────────────────────────

// True iff every face of the tet has even edge-count sum.
inline bool is_even_sum(const EdgeInts& e) {
    for (const auto& face : ALL_TET_TRIPLETS) {
        int i = face[0], j = face[1], k = face[2];
        if ((e.get(i,j) + e.get(i,k) + e.get(j,k)) % 2 != 0) return false;
    }
    return true;
}
