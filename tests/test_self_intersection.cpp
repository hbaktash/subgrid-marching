#include <catch2/catch_all.hpp>

#include "tet_generators.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/primal_reconstruction.h"
#include "subgrid_MT/cv_interpolation.h"
#include "common/file_IO.h"

#include <numeric>
#include <fstream>
#include <cmath>

// ── Edge-vs-triangle intersection test ───────────────────────────────────────

static Vector3 vcross(const Vector3& a, const Vector3& b) {
    return Vector3{a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

static double vdot(const Vector3& a, const Vector3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Möller–Trumbore: does segment (p0,p1) intersect triangle (v0,v1,v2)?
static bool segment_intersects_triangle(
    const Vector3& p0, const Vector3& p1,
    const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
    const double EPS = 1e-14;
    Vector3 dir = p1 - p0;
    Vector3 e1 = v1 - v0, e2 = v2 - v0;
    Vector3 h = vcross(dir, e2);
    double a = vdot(e1, h);
    if (a > -EPS && a < EPS) return false; // parallel
    double f = 1.0 / a;
    Vector3 s = p0 - v0;
    double u = f * vdot(s, h);
    if (u < -EPS || u > 1.0 + EPS) return false;
    Vector3 q = vcross(s, e1);
    double v = f * vdot(dir, q);
    if (v < -EPS || u + v > 1.0 + EPS) return false;
    double t = f * vdot(e2, q);
    return t > EPS && t < 1.0 - EPS;
}

// Two triangles intersect if any edge of one pierces the other
static bool triangles_intersect(
    const Vector3& a0, const Vector3& a1, const Vector3& a2,
    const Vector3& b0, const Vector3& b1, const Vector3& b2)
{
    // Edges of A vs triangle B
    if (segment_intersects_triangle(a0, a1, b0, b1, b2)) return true;
    if (segment_intersects_triangle(a1, a2, b0, b1, b2)) return true;
    if (segment_intersects_triangle(a2, a0, b0, b1, b2)) return true;
    // Edges of B vs triangle A
    if (segment_intersects_triangle(b0, b1, a0, a1, a2)) return true;
    if (segment_intersects_triangle(b1, b2, a0, a1, a2)) return true;
    if (segment_intersects_triangle(b2, b0, a0, a1, a2)) return true;
    return false;
}

// Check if two faces share a vertex
static bool faces_adjacent(const std::vector<size_t>& f1, const std::vector<size_t>& f2) {
    for (size_t a : f1)
        for (size_t b : f2)
            if (a == b) return true;
    return false;
}

// ── helpers ──────────────────────────────────────────────────────────────────

struct CombMergeGuard {
    bool prev;
    CombMergeGuard(bool val) : prev(TriangleSoup::COMB_MERGE) { TriangleSoup::COMB_MERGE = val; }
    ~CombMergeGuard() { TriangleSoup::COMB_MERGE = prev; }
};

static std::string edge_ints_str(const EdgeInts& e) {
    return std::to_string(e.get(0,1)) + "_" + std::to_string(e.get(0,2)) + "_"
         + std::to_string(e.get(0,3)) + "_" + std::to_string(e.get(1,2)) + "_"
         + std::to_string(e.get(1,3)) + "_" + std::to_string(e.get(2,3));
}

struct IntersectionPair { int face_a, face_b; };

static void save_failure(const TriangleSoup& soup, const EdgeInts& e,
                          int trial, double bulge,
                          const std::vector<IntersectionPair>& pairs,
                          const std::array<std::vector<double>,6>& ts) {
    std::string failures_dir = std::string(TEST_DATA_DIR) + "/../failures";
    std::string base = failures_dir + "/SI_" + edge_ints_str(e)
                     + "_trial" + std::to_string(trial)
                     + "_bulge" + (bulge > 0 ? "on" : "off");
    save_polygon_soup_as_obj(base + ".obj", soup.vertices, soup.faces);

    std::ofstream f(base + ".txt");
    f << "edge_ints: " << e.get(0,1) << " " << e.get(0,2) << " " << e.get(0,3)
      << " " << e.get(1,2) << " " << e.get(1,3) << " " << e.get(2,3) << "\n";
    f << "trial: " << trial << "\n";
    f << "bulge: " << bulge << "\n";
    f << "n_intersecting_pairs: " << (int)pairs.size() << "\n";
    f << "intersecting face pairs (0-indexed):\n";
    for (const auto& p : pairs)
        f << "  " << p.face_a << " " << p.face_b << "\n";
    f << "intersection t-values per edge:\n";
    for (int ei = 0; ei < 6; ++ei) {
        f << "  edge " << ALL_TET_PAIRS[ei].first << "," << ALL_TET_PAIRS[ei].second << ": ";
        for (double t : ts[ei]) f << t << " ";
        f << "\n";
    }
}

static std::vector<IntersectionPair> find_self_intersections(const TriangleSoup& soup) {
    std::vector<IntersectionPair> pairs;
    int nf = (int)soup.faces.size();
    for (int i = 0; i < nf; ++i) {
        if (soup.faces[i].size() < 3) continue;
        const auto& fi = soup.faces[i];
        const Vector3& a0 = soup.vertices[fi[0]];
        const Vector3& a1 = soup.vertices[fi[1]];
        const Vector3& a2 = soup.vertices[fi[2]];
        for (int j = i+1; j < nf; ++j) {
            if (soup.faces[j].size() < 3) continue;
            if (faces_adjacent(fi, soup.faces[j])) continue;
            const auto& fj = soup.faces[j];
            const Vector3& b0 = soup.vertices[fj[0]];
            const Vector3& b1 = soup.vertices[fj[1]];
            const Vector3& b2 = soup.vertices[fj[2]];
            if (triangles_intersect(a0, a1, a2, b0, b1, b2))
                pairs.push_back({i, j});
        }
    }
    return pairs;
}

// ── generators ───────────────────────────────────────────────────────────────

static std::vector<EdgeInts> si_inputs() {
    auto v = enumerate_up_to_symmetry(10);
    auto r = random_edge_coords(30, 10000, 10, 42);
    v.insert(v.end(), r.begin(), r.end());
    return v;
}

// ── tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("Single-tet: no self-intersections with bulge",
          "[self_intersection][single_tet][bulge]")
{
    static const auto tet_pos = standard_tet_positions();
    static const std::array<size_t,4> idx = {0, 1, 2, 3};
    CombMergeGuard guard(true);
    std::mt19937 rng(55);
    int tested = 0;

    for (const auto& e : si_inputs()) {
        auto counts = e.to_array();
        int total = std::accumulate(counts.begin(), counts.end(), 0);
        if (total == 0) continue;

        for (int trial = 0; trial < 4; ++trial) {
            auto ts = random_isect_ts(e, rng);

            EdgeOccupations occ;
            auto [open_c, scoop_c, normal_c] = boundary_comb_curves(idx, counts, occ, false);

            TriangleSoup soup = subgrid_surface(
                tet_pos, idx, ts,
                open_c, scoop_c, occ,
                true, 1e-4, false
            );
            if (soup.faces.empty()) continue;
            ++tested;

            auto si_pairs = find_self_intersections(soup);
            int n_si = (int)si_pairs.size();

            INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
                 <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3)
                 << "  trial: " << trial << "  intersections: " << n_si);

            if (n_si > 0)
                save_failure(soup, e, trial, 0.001, si_pairs, ts);

            CHECK(n_si == 0);
        }
    }
    REQUIRE(tested > 0);
}

TEST_CASE("Single-tet: no self-intersections without bulge",
          "[self_intersection][single_tet][no_bulge]")
{
    static const auto tet_pos = standard_tet_positions();
    static const std::array<size_t,4> idx = {0, 1, 2, 3};
    CombMergeGuard guard(true);
    std::mt19937 rng(55);
    int tested = 0;

    for (const auto& e : si_inputs()) {
        auto counts = e.to_array();
        int total = std::accumulate(counts.begin(), counts.end(), 0);
        if (total == 0) continue;

        for (int trial = 0; trial < 4; ++trial) {
            auto ts = random_isect_ts(e, rng);

            EdgeOccupations occ;
            auto [open_c, scoop_c, normal_c] = boundary_comb_curves(idx, counts, occ, false);

            TriangleSoup soup = subgrid_surface(
                tet_pos, idx, ts,
                open_c, scoop_c, occ,
                false, 0.0, false
            );
            if (soup.faces.empty()) continue;
            ++tested;

            auto si_pairs = find_self_intersections(soup);
            int n_si = (int)si_pairs.size();

            INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
                 <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3)
                 << "  trial: " << trial << "  intersections: " << n_si);

            if (n_si > 0)
                save_failure(soup, e, trial, 0.0, si_pairs, ts);

            CHECK(n_si == 0);
        }
    }
    REQUIRE(tested > 0);
}
