#include <catch2/catch_all.hpp>

#include "tet_generators.h"
#include "nc/NC_solver.h"

// ── helpers ───────────────────────────────────────────────────────────────────

// Combined generators: exhaustive up to max_k=3 plus random samples.
static std::vector<EdgeInts> all_inputs() {
    auto v = enumerate_up_to_symmetry(10);
    auto r = random_edge_coords(30, 10000, 10);
    v.insert(v.end(), r.begin(), r.end());
    return v;
}

// For every face {i,j,k} and every edge {a,b,c=other} of that face, check:
//   scc[a,b,c] + scc[b,a,c] + 2*scoop[c,a,b] == e[{a,b}]
// This holds universally: when triangle_ineq holds, scoop=0 so it reduces to
// the corner-coverage identity; when it fails, scoop accounts for the remainder.
static void check_edge_budget(const NormalCoordinates& nc, const EdgeInts& e) {
    for (const auto& face : ALL_TET_TRIPLETS) {
        int verts[3] = { face[0], face[1], face[2] };
        for (int ei = 0; ei < 3; ++ei) {
            int a = verts[ei], b = verts[(ei+1)%3], c = verts[(ei+2)%3];
            int lhs = nc.scc.get(a,b,c) + nc.scc.get(b,a,c)
                    + 2 * nc.scoop.get(c,a,b);
            INFO("face (" << face[0] << "," << face[1] << "," << face[2]
                 << ")  edge {" << a << "," << b << "}  lhs=" << lhs
                 << "  e=" << e.get(a,b));
            CHECK(lhs == e.get(a,b));
        }
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("NC solver: triangle-ineq holds — zero scoops and corner coverage",
          "[nc_solver][even_sum][triangle_ineq]")
{
    int tested = 0;
    for (const auto& e : all_inputs()) {
        if (!is_even_sum(e)) continue;
        NormalCoordinates nc(e);
        if (!nc.triangle_ineq) continue;
        ++tested;

        INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
             <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3));

        for (const auto& face : ALL_TET_TRIPLETS) {
            int i = face[0], j = face[1], k = face[2];
            CHECK(nc.scoop.get(i,j,k) == 0);
            CHECK(nc.scoop.get(j,i,k) == 0);
            CHECK(nc.scoop.get(k,i,j) == 0);
        }
        check_edge_budget(nc, e);
    }
    REQUIRE(tested > 0);  // guard against generator accidentally producing nothing
}

TEST_CASE("NC solver: triangle-ineq fails — scoop present and edge budget holds",
          "[nc_solver][even_sum][scoop]")
{
    int tested = 0;
    for (const auto& e : all_inputs()) {
        if (!is_even_sum(e)) continue;
        NormalCoordinates nc(e);
        if (nc.triangle_ineq) continue;
        ++tested;

        INFO("EdgeInts: " << e.get(0,1) <<" "<< e.get(0,2) <<" "<< e.get(0,3)
             <<" "<< e.get(1,2) <<" "<< e.get(1,3) <<" "<< e.get(2,3));

        bool any_scoop = false;
        for (const auto& face : ALL_TET_TRIPLETS) {
            int i = face[0], j = face[1], k = face[2];
            if (nc.scoop.get(i,j,k) + nc.scoop.get(j,i,k) + nc.scoop.get(k,i,j) > 0)
                any_scoop = true;
        }
        CHECK(any_scoop);
        check_edge_budget(nc, e);
    }
    REQUIRE(tested > 0);
}
