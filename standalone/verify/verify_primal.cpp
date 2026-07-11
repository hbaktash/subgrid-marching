// Verification harness: compares the standalone primal port against the real
// src/subgrid_MT core, in-process, over many edge-int configurations.
//
// Ground truth = the repository core. The HTML port is NOT authoritative.
// Built via the root CMake target `verify_port` (links normal_coordinates).

#include "tet_generators.h"                 // real core: EdgeInts, Vector3, generators
#include "subgrid_MT/boundary_curve.h"      // real core: boundary_comb_curves
#include "subgrid_MT/primal_reconstruction.h"  // real core: subgrid_surface
#include "common/triangle_soup.h"

#include "subgrid_mt.hpp"                    // the standalone port (namespace subgrid_mt)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <vector>

namespace {

// Match the shipped defaults (subgrid_pipeline.h): scoop mid-vertices on, bulge
// 1e-4. This exercises the simplicial-embedding Steiner path. (Corner recursion
// is now a fixed construction policy in the core, not a flag.)
constexpr bool   SCOOP_MID = true;
constexpr double BULGE     = 1e-4;

struct RefSoup { std::vector<Vector3> V; std::vector<std::vector<size_t>> F; std::vector<CombVertexKey> S; };

bool compute_reference(const std::array<Vector3,4>& pos, const std::array<size_t,4>& idx,
                       const std::array<std::vector<double>,6>& ts, RefSoup& out) {
    TriangleSoup::COMB_MERGE = true;
    auto counts = EdgeInts(ts).to_array();
    EdgeOccupations occ;
    auto [open_c, scoop_c, normal_c] = boundary_comb_curves(idx, counts, occ, /*populate_normal=*/false);
    TriangleSoup soup = subgrid_surface(pos, idx, ts, open_c, scoop_c, occ, SCOOP_MID, BULGE, /*verbose=*/false);
    out.V = soup.vertices; out.F = soup.faces; out.S = soup.signatures;
    return !open_c.empty();  // non_even
}

bool compute_port(const std::array<Vector3,4>& pos, const std::array<size_t,4>& idx,
                  const std::array<std::vector<double>,6>& ts, subgrid_mt::PrimalResult& out) {
    std::array<subgrid_mt::Vec3,4> ppos;
    for (int i = 0; i < 4; i++) ppos[i] = {pos[i].x, pos[i].y, pos[i].z};
    out = subgrid_mt::subgrid_primal(ppos, idx, ts, SCOOP_MID, BULGE);
    return out.non_even;
}

}  // namespace

int main(int argc, char** argv) {
    // Optional: dump (input, reference soup) fixtures for the Python port to verify against.
    std::ofstream fx;
    if (argc > 1) { fx.open(argv[1]); fx.precision(17); }

    auto configs = enumerate_up_to_symmetry(6);
    auto rnd = random_edge_coords(12, 3000, 6, /*seed=*/7);
    configs.insert(configs.end(), rnd.begin(), rnd.end());

    std::array<size_t,4> idx = {10, 3, 27, 6};  // non-trivial global indices to exercise remap/flip
    auto pos = standard_tet_positions();
    std::mt19937 rng(12345);

    long total = 0, compared = 0, mismatches = 0, both_threw = 0;
    double max_pos_diff = 0.0;

    for (const auto& e : configs) {
        auto counts = e.to_array();
        int sum = 0; for (int c : counts) sum += c;
        if (sum == 0) continue;
        auto ts = random_isect_ts(e, rng);
        ++total;

        RefSoup ref; subgrid_mt::PrimalResult port;
        bool ref_threw = false, port_threw = false;
        bool ref_noneven = false, port_noneven = false;
        try { ref_noneven = compute_reference(pos, idx, ts, ref); } catch (...) { ref_threw = true; }
        try { port_noneven = compute_port(pos, idx, ts, port); } catch (...) { port_threw = true; }

        if (ref_threw != port_threw) {
            printf("MISMATCH (throw): ref_threw=%d port_threw=%d  edges=[%d %d %d %d %d %d]\n",
                   ref_threw, port_threw, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);
            ++mismatches; continue;
        }
        if (ref_threw) { ++both_threw; continue; }
        ++compared;

        if (fx.is_open()) {
            fx << "C " << counts[0] << " " << counts[1] << " " << counts[2] << " "
               << counts[3] << " " << counts[4] << " " << counts[5] << "\n";
            fx << "I " << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << "\n";
            for (int e = 0; e < 6; e++) {
                fx << "T " << ts[e].size();
                for (double t : ts[e]) fx << " " << std::scientific << t;
                fx << "\n";
            }
            fx << "N " << (ref_noneven ? 1 : 0) << "\n";
            fx << "V " << ref.V.size() << "\n";
            for (const auto& v : ref.V) fx << std::scientific << v.x << " " << v.y << " " << v.z << "\n";
            fx << "F " << ref.F.size() << "\n";
            for (const auto& f : ref.F) { fx << f.size(); for (size_t id : f) fx << " " << id; fx << "\n"; }
            fx << "S " << ref.S.size() << "\n";
            for (const auto& s : ref.S)
                fx << s.i << " " << s.j << " " << s.order << " " << s.edge_int << " " << s.k
                   << " " << static_cast<int>(s.type) << "\n";
        }

        bool ok = true;
        if (ref_noneven != port_noneven) ok = false;
        if (ref.V.size() != port.vertices.size()) ok = false;
        if (ref.F.size() != port.faces.size()) ok = false;
        if (ref.S.size() != port.signatures.size()) ok = false;
        if (ok) {
            for (size_t v = 0; v < ref.V.size(); v++) {
                max_pos_diff = std::max({max_pos_diff, std::abs(ref.V[v].x - port.vertices[v].x),
                                         std::abs(ref.V[v].y - port.vertices[v].y),
                                         std::abs(ref.V[v].z - port.vertices[v].z)});
            }
            for (size_t f = 0; f < ref.F.size() && ok; f++)
                if (ref.F[f] != port.faces[f]) ok = false;
            for (size_t s = 0; s < ref.S.size() && ok; s++) {
                const auto& a = ref.S[s]; const auto& b = port.signatures[s];
                if (a.i != b.i || a.j != b.j || a.order != b.order || a.k != b.k ||
                    static_cast<int>(a.type) != static_cast<int>(b.type) || a.edge_int != b.edge_int)
                    ok = false;
            }
        }
        if (!ok) {
            ++mismatches;
            if (mismatches <= 20)
                printf("MISMATCH: edges=[%d %d %d %d %d %d] refV=%zu portV=%zu refF=%zu portF=%zu ne(r/p)=%d/%d\n",
                       counts[0], counts[1], counts[2], counts[3], counts[4], counts[5],
                       ref.V.size(), port.vertices.size(), ref.F.size(), port.faces.size(),
                       (int)ref_noneven, (int)port_noneven);
        }
    }

    printf("\n=== primal verification ===\n");
    printf("configs run:        %ld\n", total);
    printf("both threw (skip):  %ld\n", both_threw);
    printf("compared:           %ld\n", compared);
    printf("mismatches:         %ld\n", mismatches);
    printf("max position diff:  %.3e\n", max_pos_diff);
    printf("%s\n", mismatches == 0 ? "PASS" : "FAIL");
    return mismatches == 0 ? 0 : 1;
}
