#include <catch2/catch_all.hpp>

#include "subgrid_pipeline.h"
#include "dual_subgrid_pipeline.h"
#include "query/input_query_handler.h"
#include "query/cgal_queries.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"

#include <cmath>
#include <filesystem>
#include <random>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// The query cache and (later) the threaded driver only reorganize *when* work
// happens, never what it computes, so their output must match the direct path.
//
// "Match" is exact for everything topological -- vertex count, face index lists,
// and the per-tet counters -- because a mis-flipped edge or a stale cache entry
// shows up there immediately. Vertex *positions* get a float-epsilon tolerance
// for one specific, measured reason: the cache queries each grid edge once in a
// canonical direction, and neither handler is direction-symmetric (see the first
// test). Same-configuration comparisons still demand bit-for-bit equality.

namespace {

struct ParityInput {
    std::string name;
    std::function<std::unique_ptr<InputQueryHandler>()> make_handler;  // null for npz
    std::string npz_path;
    size_t res = 16;
    double pos_tol = 0.0;   // how far the cache may move a vertex, per handler
    // Whether the handler returns the same crossing *count* in both directions.
    // False for SDF: its march can straddle a pair of nearby crossings (see the
    // direction test), so canonicalizing the query direction can change topology.
    bool direction_stable = true;
    bool is_npz() const { return !npz_path.empty(); }
};

std::vector<std::string> files_with_ext(const std::string& dir, const std::vector<std::string>& exts) {
    std::vector<std::string> files;
    if (!fs::exists(dir)) return files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        auto ext = entry.path().extension().string();
        for (const auto& e : exts)
            if (ext == e) { files.push_back(entry.path().string()); break; }
    }
    std::sort(files.begin(), files.end());
    return files;
}

#ifdef HAVE_CGAL
// Exact EPECK queries. Only reachable in a -DSUBGRID_WITH_CGAL=ON build, and the
// reason the parity/thread tests cover it is that its lazy kernel reference-counts
// its handles: those counters are non-atomic unless CGAL_HAS_THREADS is defined,
// so this is the input that would expose a missing define under ThreadSanitizer.
std::unique_ptr<InputQueryHandler> cgal_handler(const std::string& path) {
    geometrycentral::surface::SimplePolygonMesh simple_mesh;
    simple_mesh.readMeshFromFile(path, "");
    auto pre = preprocess_input_mesh(simple_mesh);
    return std::make_unique<CGALQueryHandler>(pre.positions, pre.polygons);
}
#endif

std::unique_ptr<InputQueryHandler> mesh_handler(const std::string& path) {
    geometrycentral::surface::SimplePolygonMesh simple_mesh;
    simple_mesh.readMeshFromFile(path, "");
    auto pre = preprocess_input_mesh(simple_mesh);
    return std::make_unique<MeshQueryHandler>(pre.positions, pre.polygons);
}

// Measured worst-case direction asymmetry of a single edge query, over the
// 400-segment sweep in the direction test below:
//   mesh (FCPW, float internals) ......... 3.3e-6
//   SDF  (marching, step 1e-2) ........... 2.1e-4   <- set by the march step
//   precomputed (npz) .................... exact, it already canonicalizes
// These bound how far the cache can move a vertex relative to the uncached path.
constexpr double MESH_DIRECTION_TOL = 1e-5;
constexpr double SDF_DIRECTION_TOL  = 1e-3;
// CGAL/EPECK is direction-symmetric in exact arithmetic -- the rational t and
// 1-t are exact. The whole residue is the final rounding to double: to_double(t)
// and 1.0 - to_double(1-t) round independently and can land one ULP apart
// (measured 1.11e-16 on t), which the interpolation and centroid arithmetic in
// the construction then carries into position (measured 1.55e-15, ~7 ULPs on a
// coordinate in [-1,1]). Slack for that, and nothing that could hide a real
// disagreement -- the next-least-exact handler is 1e4 times coarser.
constexpr double CGAL_DIRECTION_TOL = 1e-13;

// One representative of each input family; the cache is orthogonal to the
// construction, so breadth over resolutions matters more than over models.
std::vector<ParityInput> parity_inputs() {
    std::vector<ParityInput> inputs;
    for (size_t res : {8u, 16u, 24u}) {
        inputs.push_back({"sdf:Sphere", [] { return std::make_unique<SDFQueryHandler>("Sphere", 1e-2f); }, "", res, SDF_DIRECTION_TOL, false});
        inputs.push_back({"sdf:Cables", [] { return std::make_unique<SDFQueryHandler>("Cables", 1e-2f); }, "", res, SDF_DIRECTION_TOL, false});
        auto meshes = files_with_ext(TEST_DATA_DIR, {".obj", ".ply", ".stl", ".off"});
        if (!meshes.empty()) {
            const auto path = meshes.front();
            inputs.push_back({"mesh:" + fs::path(path).stem().string(),
                              [path] { return mesh_handler(path); }, "", res, MESH_DIRECTION_TOL, true});
#ifdef HAVE_CGAL
            // Exact predicates *and* exact constructions, so unlike FCPW this one
            // is direction-symmetric up to the final rounding to double.
            inputs.push_back({"cgal:" + fs::path(path).stem().string(),
                              [path] { return cgal_handler(path); }, "", res, CGAL_DIRECTION_TOL, true});
#endif
        }
    }
    for (const auto& path : files_with_ext(TEST_NPZ_DIR, {".npz"}))
        inputs.push_back({"npz:" + fs::path(path).stem().string(), nullptr, path, 0, 0.0, true});
    return inputs;
}

struct CombMergeGuard {
    bool prev;
    CombMergeGuard(bool val) : prev(TriangleSoup::COMB_MERGE) { TriangleSoup::COMB_MERGE = val; }
    ~CombMergeGuard() { TriangleSoup::COMB_MERGE = prev; }
};

// Topology is compared exactly -- same vertex count, same face index lists, in
// the same order. Positions take a tolerance only because the cache canonicalizes
// each edge's query direction and neither handler is direction-symmetric (see the
// direction test below); `pos_tol = 0` demands bit-for-bit equality and is what
// every same-configuration comparison uses.
void require_same_soup(const TriangleSoup& a, const TriangleSoup& b, double pos_tol = 0.0) {
    REQUIRE(a.vertices.size() == b.vertices.size());
    REQUIRE(a.faces.size() == b.faces.size());
    for (size_t f = 0; f < a.faces.size(); ++f)
        REQUIRE(a.faces[f] == b.faces[f]);
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        const Vector3 d = a.vertices[i] - b.vertices[i];
        INFO("vertex " << i << " differs by " << d.norm());
        REQUIRE(std::abs(d.x) <= pos_tol);
        REQUIRE(std::abs(d.y) <= pos_tol);
        REQUIRE(std::abs(d.z) <= pos_tol);
    }
}

SubgridPipelineResult run_primal(const ParityInput& in, bool cache, int threads = 1,
                                bool canonical = false) {
    SubgridPipelineOpts opts;
    opts.query_cache = cache;
    opts.num_threads = threads;
    opts.canonical_queries = canonical;
    if (in.is_npz()) return run_subgrid_pipeline_npz(in.npz_path, opts);
    auto handler = in.make_handler();
    return run_subgrid_pipeline(*handler, in.res, opts);
}

DualSubgridPipelineResult run_dual(const ParityInput& in, bool cache, int threads = 1,
                                  bool canonical = false) {
    DualSubgridPipelineOpts opts;
    opts.query_cache = cache;
    opts.num_threads = threads;
    opts.canonical_queries = canonical;
    opts.no_normal = in.is_npz();  // the bundled npz files may carry no normals
    if (in.is_npz()) return run_dual_subgrid_pipeline_npz(in.npz_path, opts);
    auto handler = in.make_handler();
    return run_dual_subgrid_pipeline(*handler, in.res, opts);
}

}  // namespace


TEST_CASE("Handlers: edge queries are direction-asymmetric at float precision",
          "[parity][query_cache][direction]") {
    // Neither handler is direction-symmetric: FCPW works in float internally, and
    // the SDF march samples at different points depending on which end it starts
    // from. So querying b->a does NOT give exactly 1 - t of querying a->b.
    //
    // This is pre-existing, not something the cache introduced: without a cache
    // each tet queries its edges in its own local direction, so two tets sharing
    // a grid edge already disagree about the crossing position at this scale.
    // The cache queries once, canonically, which makes them agree -- and shifts
    // positions by up to this much relative to the uncached path.
    //
    // What must hold is that the disagreement stays at float epsilon and that the
    // crossing *count* is direction-independent (counts drive all the topology).
    // Swept over many random segments -- a single probe badly understates the
    // SDF march, whose bracketing interval shifts by up to a full step when the
    // direction flips.
    auto probe = [&](InputQueryHandler& h, const char* name,
                     double max_t_asymmetry, double max_count_mismatch_rate) {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-0.6, 0.6);

        double worst = 0.0;
        size_t crossings = 0, count_mismatches = 0;
        for (int trial = 0; trial < 400; ++trial) {
            const Vector3 a{u(rng), u(rng), u(rng)};
            const Vector3 b{u(rng), u(rng), u(rng)};

            std::vector<double> fwd_ts, rev_ts;
            std::vector<Vector3> fwd_n, rev_n;
            h.query_edge(0, 1, a, b, fwd_ts, fwd_n, false, true);
            h.query_edge(1, 0, b, a, rev_ts, rev_n, false, true);

            if (fwd_ts.size() != rev_ts.size()) { count_mismatches++; continue; }
            crossings += fwd_ts.size();
            for (size_t k = 0; k < fwd_ts.size(); ++k)
                worst = std::max(worst, std::abs(fwd_ts[k] - (1.0 - rev_ts[rev_ts.size() - 1 - k])));
        }

        const double mismatch_rate = double(count_mismatches) / 400.0;
        INFO("handler=" << name << "  worst |t_fwd - (1 - t_rev)| = " << worst
             << "  crossings=" << crossings
             << "  count mismatches=" << count_mismatches << " (" << mismatch_rate * 100 << "%)");
        REQUIRE(crossings > 0);   // the sweep must actually hit the surface
        // Crossing counts drive every topological decision, so where they are
        // direction-stable the cache cannot change topology at all.
        CHECK(mismatch_rate <= max_count_mismatch_rate);
        CHECK(worst <= max_t_asymmetry);
    };

    SECTION("SDF") {
        // Root-finding by marching: the sample grid shifts with direction, so the
        // linear interpolation lands on a slightly different root. Much larger
        // than pure float noise, and it scales with the march step.
        SDFQueryHandler h("Sphere", 1e-2f);
        probe(h, "SDFQueryHandler", SDF_DIRECTION_TOL, /*counts stable=*/0.0);
    }

    SECTION("SDF, high-genus") {
        // Cables has many thin features, so a step can straddle a pair of nearby
        // crossings -- the case where direction can change the crossing *count*,
        // not just its position.
        SDFQueryHandler h("Cables", 1e-2f);
        probe(h, "SDFQueryHandler/Cables", 1e-2, /*measured 2.5% disagree=*/0.05);
    }

#ifdef HAVE_CGAL
    SECTION("CGAL / EPECK") {
        // Exact predicates AND exact constructions: t and 1-t are exact rationals,
        // so the only residue is rounding each to double independently at the end.
        auto meshes = files_with_ext(TEST_DATA_DIR, {".obj", ".ply", ".stl", ".off"});
        if (meshes.empty()) { WARN("no mesh inputs; set TEST_DATA_DIR"); return; }
        auto h = cgal_handler(meshes.front());
        probe(*h, "CGALQueryHandler", CGAL_DIRECTION_TOL, /*counts stable=*/0.0);
    }
#endif

    SECTION("mesh") {
        // Exact ray/triangle hit, so the only asymmetry is FCPW's internal float.
        auto meshes = files_with_ext(TEST_DATA_DIR, {".obj", ".ply", ".stl", ".off"});
        if (meshes.empty()) { WARN("no mesh inputs; set TEST_DATA_DIR"); return; }
        auto h = mesh_handler(meshes.front());
        probe(*h, "MeshQueryHandler", MESH_DIRECTION_TOL, /*counts stable=*/0.0);
    }
}


TEST_CASE("Query cache: primal output matches the direct path",
          "[parity][query_cache]") {
    CombMergeGuard guard(true);
    for (const auto& in : parity_inputs()) {
        if (!in.direction_stable) continue;   // covered by the SDF test below
        const auto baseline = run_primal(in, false);

        for (bool cache : {true}) {
            INFO("input=" << in.name << " res=" << in.res
                 << " cache=" << cache);
            const auto cached = run_primal(in, cache);
            REQUIRE(cached.non_zero_tets   == baseline.non_zero_tets);
            REQUIRE(cached.non_even_tets   == baseline.non_even_tets);
            REQUIRE(cached.non_normal_tets == baseline.non_normal_tets);
            REQUIRE(cached.total_tets      == baseline.total_tets);
            require_same_soup(cached.soup, baseline.soup, in.pos_tol);
        }
    }
}


TEST_CASE("Query cache: dual output matches the direct path",
          "[parity][query_cache]") {
    for (const auto& in : parity_inputs()) {
        if (!in.direction_stable) continue;   // covered by the SDF test below
        const auto baseline = run_dual(in, false);

        for (bool cache : {true}) {
            INFO("input=" << in.name << " res=" << in.res
                 << " cache=" << cache);
            const auto cached = run_dual(in, cache);
            REQUIRE(cached.non_zero_tets   == baseline.non_zero_tets);
            REQUIRE(cached.non_even_tets   == baseline.non_even_tets);
            REQUIRE(cached.non_normal_tets == baseline.non_normal_tets);
            require_same_soup(cached.soup, baseline.soup, in.pos_tol);
            REQUIRE(cached.soup.dual_positions.size() == baseline.soup.dual_positions.size());
            for (size_t k = 0; k < cached.soup.dual_positions.size(); ++k) {
                const Vector3 d = cached.soup.dual_positions[k] - baseline.soup.dual_positions[k];
                INFO("dual point " << k << " differs by " << d.norm());
                // The QEF amplifies the input shift, so allow a looser bound here
                // than on the boundary vertices that feed it.
                REQUIRE(d.norm() <= 100.0 * in.pos_tol + 1e-12);
            }
            REQUIRE(cached.soup.faces_per_edge == baseline.soup.faces_per_edge);
        }
    }
}


TEST_CASE("Query cache: SDF input shifts slightly, because its march is direction-unstable",
          "[parity][query_cache][sdf]") {
    // The SDF handler's crossing count is not direction-independent: on thin
    // features a march step can straddle a pair of nearby crossings, and the
    // direction test measures ~2.5% of segments disagreeing on `Cables`.
    //
    // Without a cache each tet queries a shared grid edge in its own direction,
    // so two tets can disagree about how many times the surface crosses the edge
    // they share. The cache queries once, so they always agree -- which is more
    // self-consistent, but not identical to the uncached result. Hence no exact
    // topological claim here; what we do require is that the difference stays
    // marginal and that the cache never makes the even-sum condition worse.
    CombMergeGuard guard(true);
    for (const auto& in : parity_inputs()) {
        if (in.direction_stable) continue;

        const auto baseline = run_primal(in, false);
        for (bool cache : {true}) {
            INFO("input=" << in.name << " res=" << in.res
                 << " cache=" << cache);
            const auto cached = run_primal(in, cache);

            REQUIRE(cached.total_tets == baseline.total_tets);
            // Consistent shared-edge counts cannot introduce even-sum violations.
            CHECK(cached.non_even_tets <= baseline.non_even_tets);

            auto near = [](size_t a, size_t b) {
                const double lo = double(std::min(a, b)), hi = double(std::max(a, b));
                return hi <= lo * 1.01 + 8;   // within 1% (plus slack at tiny sizes)
            };
            INFO("verts " << cached.soup.vertices.size() << " vs " << baseline.soup.vertices.size()
                 << ", faces " << cached.soup.faces.size() << " vs " << baseline.soup.faces.size());
            CHECK(near(cached.soup.vertices.size(), baseline.soup.vertices.size()));
            CHECK(near(cached.soup.faces.size(), baseline.soup.faces.size()));
        }
    }
}

TEST_CASE("Canonical queries: the cache becomes bit-exact against the direct path",
          "[parity][canonical]") {
    // The cache always queries an edge once, in the canonical min->max direction.
    // The uncached path normally queries it once per incident tet, in each tet's
    // own direction -- and no handler is exactly direction-symmetric, which is the
    // entire reason the cache tests above need a position tolerance and why SDF
    // input needs a separate, weaker test.
    //
    // canonical_queries makes the uncached path ask the same way. Then the two
    // must agree exactly, for every input including SDF, with no tolerance and no
    // topological wiggle room. It costs no extra queries -- only a reversal on the
    // roughly half of tet edges that run high-index to low.
    CombMergeGuard guard(true);
    for (const auto& in : parity_inputs()) {
        const auto baseline = run_primal(in, false, 1, /*canonical=*/true);

        for (bool cache : {true}) {
            INFO("primal input=" << in.name << " res=" << in.res
                 << " cache=" << cache);
            const auto cached = run_primal(in, cache, 1, /*canonical=*/true);
            REQUIRE(cached.non_zero_tets   == baseline.non_zero_tets);
            REQUIRE(cached.non_even_tets   == baseline.non_even_tets);
            REQUIRE(cached.non_normal_tets == baseline.non_normal_tets);
            require_same_soup(cached.soup, baseline.soup, 0.0);
        }
    }

    for (const auto& in : parity_inputs()) {
        const auto baseline = run_dual(in, false, 1, /*canonical=*/true);
        for (bool cache : {true}) {
            INFO("dual input=" << in.name << " res=" << in.res
                 << " cache=" << cache);
            const auto cached = run_dual(in, cache, 1, /*canonical=*/true);
            require_same_soup(cached.soup, baseline.soup, 0.0);
            REQUIRE(cached.soup.faces_per_edge == baseline.soup.faces_per_edge);
            for (size_t k = 0; k < cached.soup.dual_positions.size(); ++k)
                REQUIRE(cached.soup.dual_positions[k] == baseline.soup.dual_positions[k]);
        }
    }
}


TEST_CASE("Threads: primal output is bit-identical regardless of thread count",
          "[parity][threads]") {
    // The strongest guarantee in this branch, and the reason the merge stays
    // sequential and in tet order: threading changes only *when* work happens.
    // Unlike the cache, it does not touch query direction, so there is nothing
    // to excuse here -- exact equality or the driver is wrong.
    CombMergeGuard guard(true);
    for (const auto& in : parity_inputs()) {
        for (bool cache : {false, true}) {
            const auto serial = run_primal(in, cache, 1);
            for (int threads : {2, 4, 8}) {
                INFO("input=" << in.name << " res=" << in.res << " threads=" << threads
                     << " cache=" << cache);
                const auto par = run_primal(in, cache, threads);
                REQUIRE(par.non_zero_tets   == serial.non_zero_tets);
                REQUIRE(par.non_even_tets   == serial.non_even_tets);
                REQUIRE(par.non_normal_tets == serial.non_normal_tets);
                REQUIRE(par.total_tets      == serial.total_tets);
                require_same_soup(par.soup, serial.soup, 0.0);
            }
        }
    }
}


TEST_CASE("Threads: dual output is bit-identical regardless of thread count",
          "[parity][threads]") {
    for (const auto& in : parity_inputs()) {
        for (bool cache : {false, true}) {
            const auto serial = run_dual(in, cache, 1);
            for (int threads : {2, 8}) {
                INFO("input=" << in.name << " res=" << in.res << " threads=" << threads
                     << " cache=" << cache);
                const auto par = run_dual(in, cache, threads);
                REQUIRE(par.non_zero_tets   == serial.non_zero_tets);
                REQUIRE(par.non_even_tets   == serial.non_even_tets);
                REQUIRE(par.non_normal_tets == serial.non_normal_tets);
                require_same_soup(par.soup, serial.soup, 0.0);
                REQUIRE(par.soup.faces_per_edge == serial.soup.faces_per_edge);
                REQUIRE(par.soup.dual_positions.size() == serial.soup.dual_positions.size());
                for (size_t k = 0; k < par.soup.dual_positions.size(); ++k)
                    REQUIRE(par.soup.dual_positions[k] == serial.soup.dual_positions[k]);
            }
        }
    }
}


TEST_CASE("Query cache: slab issues one query per unique grid edge",
          "[parity][query_cache]") {
    // The whole point of the cache. Without it the pipeline issues 6 queries per
    // tet = 30n^3; the grid has only 3n(n+1)^2 + 3n^2(n+1) distinct edges.
    for (size_t n : {4u, 8u, 12u}) {
        TetGridRange grid(n, true);
        SDFQueryHandler handler("Sphere", 1e-2f);

      {
        auto cache = make_edge_cache(true, grid);
        REQUIRE(cache != nullptr);

        std::array<std::vector<double>, 6> ts;
        std::array<std::vector<Vector3>, 6> normals;
        const size_t per_slab = grid.tets_per_slab();
        for (size_t flat = 0; flat < grid.total_tet_count(); ++flat) {
            const size_t chunk = flat / per_slab;
            if (flat % per_slab == 0)
                cache->begin_chunk(handler, chunk, false, false);
            auto tet = grid.tet_at(flat);
            cache->query_tet(handler, tet.indices, tet.positions, ts, normals, false, false);
        }

        const size_t unique_edges = 3 * n * (n + 1) * (n + 1) + 3 * n * n * (n + 1);
        INFO("n=" << n);
        REQUIRE(cache->edge_queries() == unique_edges);
        // ...against what the uncached path would have issued.
        REQUIRE(6 * grid.total_tet_count() == 30 * n * n * n);
      }
    }
}
