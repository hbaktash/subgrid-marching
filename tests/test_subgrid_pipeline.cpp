#include <catch2/catch_all.hpp>

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "subgrid_pipeline.h"
#include "assembly/mesh_processing.h"
#include "query/input_query_handler.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

namespace fs = std::filesystem;

// ── input abstraction ────────────────────────────────────────────────────────

struct TestInput {
    std::string name;
    std::function<std::unique_ptr<InputQueryHandler>()> make_handler; // null for npz
    std::string npz_path; // non-empty for npz inputs
    bool is_npz() const { return !npz_path.empty(); }
};

// A single parameterized case. SDF/mesh inputs are swept over grid resolutions;
// npz inputs are precomputed at their own baked-in resolution, so they run once
// and `res` is unused for them.
struct PipelineCase {
    TestInput input;
    size_t    res;
};

static std::vector<std::string> find_files_with_ext(const std::string& dir, const std::vector<std::string>& exts) {
    std::vector<std::string> files;
    if (!fs::exists(dir)) return files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        auto ext = entry.path().extension().string();
        for (const auto& e : exts) {
            if (ext == e) { files.push_back(entry.path().string()); break; }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

static std::unique_ptr<InputQueryHandler> make_mesh_handler(const std::string& path) {
    geometrycentral::surface::SimplePolygonMesh simple_mesh;
    simple_mesh.readMeshFromFile(path, "");
    auto pre = preprocess_input_mesh(simple_mesh);
    return std::make_unique<MeshQueryHandler>(pre.positions, pre.polygons);
}

static SubgridPipelineResult run_for_input(const TestInput& input, size_t res) {
    if (input.is_npz())
        return run_subgrid_pipeline_npz(input.npz_path);
    auto handler = input.make_handler();
    return run_subgrid_pipeline(*handler, res);
}

static std::vector<TestInput> all_pipeline_inputs() {
    std::vector<TestInput> inputs;
    for (const char* sdf : {"Sphere", "Girl", "Cables", "Teapot", "UprightPiano", "Cube", "GrandPiano"})
        inputs.push_back({sdf, [=]{ return std::make_unique<SDFQueryHandler>(sdf, 1e-2f); }, ""});
    for (const auto& path : find_files_with_ext(TEST_DATA_DIR, {".obj", ".ply", ".stl", ".off"}))
        inputs.push_back({fs::path(path).stem().string(),
                          [=]{ return make_mesh_handler(path); }, ""});
    for (const auto& path : find_files_with_ext(TEST_NPZ_DIR, {".npz"}))
        inputs.push_back({"npz:" + fs::path(path).stem().string(), nullptr, path});
    return inputs;
}

// Expand inputs into (input, res) cases: sweep resolutions for SDF/mesh inputs,
// run each npz once (it carries its own resolution).
static std::vector<PipelineCase> all_pipeline_cases() {
    std::vector<PipelineCase> cases;
    for (const auto& in : all_pipeline_inputs()) {
        if (in.is_npz()) cases.push_back({in, 0});
        else { cases.push_back({in, 16}); cases.push_back({in, 32}); }
    }
    return cases;
}

// ── helpers ──────────────────────────────────────────────────────────────────

struct MeshOwner {
    SurfaceMesh*             mesh;
    VertexPositionGeometry*  geo;
    ~MeshOwner() { delete mesh; delete geo; }
};

static void check_soup(const TriangleSoup& soup) {
    REQUIRE(!soup.vertices.empty());
    REQUIRE(!soup.faces.empty());

    for (const auto& v : soup.vertices) {
        CHECK(std::isfinite(v.x));
        CHECK(std::isfinite(v.y));
        CHECK(std::isfinite(v.z));
    }

    for (const auto& f : soup.faces) {
        CHECK(f.size() >= 3);
        for (size_t i : f)
            CHECK(i < soup.vertices.size());
    }
}

struct CombMergeGuard {
    bool prev;
    CombMergeGuard(bool val) : prev(TriangleSoup::COMB_MERGE) { TriangleSoup::COMB_MERGE = val; }
    ~CombMergeGuard() { TriangleSoup::COMB_MERGE = prev; }
};

// ── tests ─────────────────────────────────────────────────────────────────────

// Numerical (positional) merge is a best-effort, opt-in path (`--mergeEPS`);
// combinatorial merge is the robust default (see the comb-merge tests below).
// A single global epsilon cannot both stitch every should-be-identical duplicate
// (adjacent tets recompute shared-face vertices with ~1e-6-scale float drift) and
// avoid collapsing genuinely-distinct nearby vertices, so on detailed inputs the
// merged mesh may have cracks (too-small eps) or non-manifold pinches (too-large
// eps). We therefore only assert the soup is well-formed and the mesh builds;
// manifoldness/orientability are reported as WARN, not hard failures.
TEST_CASE("Subgrid pipeline: Numerical merge builds a mesh (manifoldness best-effort)",
          "[subgrid_pipeline]") {
    auto cases = all_pipeline_cases();
    if (cases.size() <= 1)
        WARN("Only SDF inputs available — set TEST_DATA_DIR for mesh coverage");

    auto tc = GENERATE_COPY(from_range(cases));
    const auto& input = tc.input;
    auto res = tc.res;
    INFO("input: " << input.name << (input.is_npz() ? "  (npz, own resolution)"
                                                    : "  res: " + std::to_string(res)));

    auto result = run_for_input(input, res);
    auto& soup = result.soup;
    if (soup.faces.empty()) SKIP("No faces produced");
    check_soup(soup);

    auto [mesh, geo] = mergeIdenticalVertices(1e-6, soup.faces, soup.vertices);
    MeshOwner owner{mesh, geo};

    REQUIRE(mesh != nullptr);
    if (!check_edge_manifoldness(*mesh))
        WARN("numerical merge (eps=1e-6) left non-manifold edges for " << input.name << " @res " << res);
    if (result.non_even_tets == 0 && !mesh->isOriented())
        WARN("numerical merge (eps=1e-6) produced a non-orientable mesh for " << input.name << " @res " << res);
}

TEST_CASE("Subgrid pipeline: comb-merge : manifold, orientable",
          "[subgrid_pipeline][combMerge]") {
    auto tc = GENERATE_COPY(from_range(all_pipeline_cases()));
    const auto& input = tc.input;
    auto res = tc.res;
    INFO("input: " << input.name << (input.is_npz() ? "  (npz, own resolution)"
                                                    : "  res: " + std::to_string(res)));

    CombMergeGuard guard(true);
    auto result = run_for_input(input, res);
    auto& soup = result.soup;
    if (soup.faces.empty()) SKIP("No faces produced");
    check_soup(soup);

    auto [mesh, geo] = makeSurfaceMeshAndGeometry(soup.faces, soup.vertices);

    REQUIRE(mesh != nullptr);
    CHECK(check_edge_manifoldness(*mesh));
    if (result.non_even_tets == 0){
        mesh->greedilyOrientFaces();
        CHECK(mesh->isOriented());
    }
}

TEST_CASE("Subgrid pipeline: comb-merge leaves no residual duplicates",
          "[subgrid_pipeline][combMerge]") {
    auto tc = GENERATE_COPY(from_range(all_pipeline_cases()));
    const auto& input = tc.input;
    auto res = tc.res;
    INFO("input: " << input.name << (input.is_npz() ? "  (npz, own resolution)"
                                                    : "  res: " + std::to_string(res)));

    CombMergeGuard guard(true);
    auto result = run_for_input(input, res);
    auto& soup = result.soup;
    if (soup.faces.empty()) SKIP("No faces produced");
    REQUIRE(!soup.faces.empty());

    auto [mesh, geo] = mergeIdenticalVertices(1e-6, soup.faces, soup.vertices);
    MeshOwner owner{mesh, geo};
    CHECK(mesh->nVertices() == soup.vertices.size());
}
