// Verification harness: compares the standalone DUAL port against the real
// src/subgrid_MT/dual core, in-process, over many edge-int configurations.
// Ground truth = the repository core. Built via root CMake target `verify_port_dual`.

#include "tet_generators.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/dual/dual_construction.h"
#include "common/triangle_soup.h"

#include "subgrid_mt.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <vector>

namespace {
constexpr double REG_ALPHA = 0.1;  // matches the shipped default (dual_subgrid_pipeline.h)

// Random unit normals per intersection, shared between reference and port.
std::array<std::vector<std::array<double, 3>>, 6> random_normals(const std::array<std::vector<double>, 6>& ts,
                                                                 std::mt19937& rng) {
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    std::array<std::vector<std::array<double, 3>>, 6> nrm;
    for (int e = 0; e < 6; e++) {
        for (size_t k = 0; k < ts[e].size(); k++) {
            double x = d(rng), y = d(rng), z = d(rng);
            double n = std::sqrt(x * x + y * y + z * z);
            if (n < 1e-9) { x = 1; y = 0; z = 0; n = 1; }
            nrm[e].push_back({x / n, y / n, z / n});
        }
    }
    return nrm;
}
}  // namespace

int main(int argc, char** argv) {
    std::ofstream fx;
    if (argc > 1) { fx.open(argv[1]); fx.precision(17); }

    bool project = (argc > 2 && std::string(argv[2]) == "project");

    auto configs = enumerate_up_to_symmetry(6);
    auto rnd = random_edge_coords(12, 3000, 6, /*seed=*/7);
    configs.insert(configs.end(), rnd.begin(), rnd.end());

    std::array<size_t, 4> idx = {10, 3, 27, 6};
    auto pos = standard_tet_positions();
    std::array<subgrid_mt::Vec3, 4> ppos;
    for (int i = 0; i < 4; i++) ppos[i] = {pos[i].x, pos[i].y, pos[i].z};
    std::mt19937 rng(999);

    long total = 0, compared = 0, mismatches = 0, both_empty = 0;
    double max_pos_diff = 0.0, max_dual_diff = 0.0;

    for (const auto& e : configs) {
        auto counts = e.to_array();
        int sum = 0; for (int c : counts) sum += c;
        if (sum == 0) continue;
        auto ts = random_isect_ts(e, rng);
        auto nrm = random_normals(ts, rng);
        ++total;

        // typed normals
        std::array<std::vector<Vector3>, 6> ref_nrm;
        std::array<std::vector<subgrid_mt::Vec3>, 6> port_nrm;
        for (int ee = 0; ee < 6; ee++)
            for (auto& v : nrm[ee]) { ref_nrm[ee].push_back(Vector3{v[0], v[1], v[2]}); port_nrm[ee].push_back({v[0], v[1], v[2]}); }

        // reference
        TriangleSoup::COMB_MERGE = true;
        EdgeOccupations occ;
        auto [open_c, scoop_c, normal_c] = boundary_comb_curves(idx, counts, occ, /*populate_normal=*/true);
        TriangleSoup ref = dual_subgrid_surface(pos, idx, ts, ref_nrm, open_c, scoop_c, normal_c, REG_ALPHA, project);

        // port
        subgrid_mt::DualResult port = subgrid_mt::subgrid_dual(ppos, idx, ts, port_nrm, REG_ALPHA, project);

        if (ref.faces.empty() && port.faces.empty()) { ++both_empty; continue; }
        ++compared;

        bool ok = true;
        if (ref.vertices.size() != port.vertices.size()) ok = false;
        if (ref.faces.size() != port.faces.size()) ok = false;
        if (ref.faces_per_edge.size() != port.faces_per_edge.size()) ok = false;
        if (ref.dual_positions.size() != port.dual_positions.size()) ok = false;
        if (ref.signatures.size() != port.signatures.size()) ok = false;
        if (ok) {
            for (size_t v = 0; v < ref.vertices.size(); v++)
                max_pos_diff = std::max({max_pos_diff, std::abs(ref.vertices[v].x - port.vertices[v].x),
                                         std::abs(ref.vertices[v].y - port.vertices[v].y),
                                         std::abs(ref.vertices[v].z - port.vertices[v].z)});
            for (size_t f = 0; f < ref.faces.size() && ok; f++)
                if (ref.faces[f] != port.faces[f]) ok = false;
            for (size_t f = 0; f < ref.faces_per_edge.size() && ok; f++) {
                if (ref.faces_per_edge[f].size() != port.faces_per_edge[f].size()) { ok = false; break; }
                for (size_t s = 0; s < ref.faces_per_edge[f].size(); s++)
                    if (ref.faces_per_edge[f][s] != port.faces_per_edge[f][s]) { ok = false; break; }
            }
            for (size_t d = 0; d < ref.dual_positions.size(); d++) {
                const auto& a = ref.dual_positions[d]; const auto& b = port.dual_positions[d];
                bool ad = a.isDefined(), bd = b.isDefined();
                if (ad != bd) { ok = false; break; }
                if (ad) max_dual_diff = std::max({max_dual_diff, std::abs(a.x - b.x), std::abs(a.y - b.y), std::abs(a.z - b.z)});
            }
            for (size_t s = 0; s < ref.signatures.size() && ok; s++) {
                const auto& a = ref.signatures[s]; const auto& b = port.signatures[s];
                if (a.i != b.i || a.j != b.j || a.order != b.order || a.k != b.k ||
                    static_cast<int>(a.type) != static_cast<int>(b.type) || a.edge_int != b.edge_int)
                    ok = false;
            }
        }
        if (!ok) {
            ++mismatches;
            if (mismatches <= 20)
                printf("MISMATCH edges=[%d %d %d %d %d %d] refV=%zu portV=%zu refF=%zu portF=%zu\n",
                       counts[0], counts[1], counts[2], counts[3], counts[4], counts[5],
                       ref.vertices.size(), port.vertices.size(), ref.faces.size(), port.faces.size());
        }

        if (fx.is_open()) {
            fx << "C " << counts[0] << " " << counts[1] << " " << counts[2] << " "
               << counts[3] << " " << counts[4] << " " << counts[5] << "\n";
            fx << "I " << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << "\n";
            for (int ee = 0; ee < 6; ee++) { fx << "T " << ts[ee].size(); for (double t : ts[ee]) fx << " " << t; fx << "\n"; }
            for (int ee = 0; ee < 6; ee++) {
                fx << "M " << nrm[ee].size();
                for (auto& v : nrm[ee]) fx << " " << v[0] << " " << v[1] << " " << v[2];
                fx << "\n";
            }
            fx << "V " << ref.vertices.size() << "\n";
            for (const auto& v : ref.vertices) fx << v.x << " " << v.y << " " << v.z << "\n";
            fx << "F " << ref.faces.size() << "\n";
            for (const auto& f : ref.faces) { fx << f.size(); for (size_t id : f) fx << " " << id; fx << "\n"; }
            fx << "E " << ref.faces_per_edge.size() << "\n";
            for (const auto& fpe : ref.faces_per_edge) {
                fx << fpe.size();
                for (const auto& t : fpe) fx << " " << t[0] << " " << t[1] << " " << t[2];
                fx << "\n";
            }
            fx << "D " << ref.dual_positions.size() << "\n";
            for (const auto& d : ref.dual_positions) {
                if (d.isDefined()) fx << d.x << " " << d.y << " " << d.z << "\n";
                else fx << "nan nan nan\n";
            }
            fx << "S " << ref.signatures.size() << "\n";
            for (const auto& s : ref.signatures)
                fx << s.i << " " << s.j << " " << s.order << " " << s.edge_int << " " << s.k << " " << static_cast<int>(s.type) << "\n";
        }
    }

    printf("\n=== dual verification (project=%d) ===\n", (int)project);
    printf("configs run:        %ld\n", total);
    printf("both empty (skip):  %ld\n", both_empty);
    printf("compared:           %ld\n", compared);
    printf("mismatches:         %ld\n", mismatches);
    printf("max position diff:  %.3e\n", max_pos_diff);
    printf("max dual pos diff:  %.3e\n", max_dual_diff);
    printf("%s\n", mismatches == 0 ? "PASS" : "FAIL");
    return mismatches == 0 ? 0 : 1;
}
