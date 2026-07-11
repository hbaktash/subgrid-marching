#include <catch2/catch_all.hpp>

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "dual_subgrid_pipeline.h"
#include "assembly/dual_mesh_assembly.h"
#include "query/input_query_handler.h"
#include "common/npz_reader.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"

namespace fs = std::filesystem;

// ── input abstraction ────────────────────────────────────────────────────────

struct TestInput {
    std::string name;
    std::function<std::unique_ptr<InputQueryHandler>()> make_handler; // null for npz
    std::string npz_path;
    bool is_npz() const { return !npz_path.empty(); }
};

// A single parameterized case. SDF/mesh inputs are swept over grid resolutions;
// npz inputs are precomputed at their own baked-in resolution, so they run once
// and `res` is unused for them.
struct DualCase {
    TestInput input;
    size_t    res;
};

// The dual QEF solve needs a surface normal per intersection; npz files without
// an `isect_normals` array cannot drive the dual pipeline.
static bool npz_has_normals(const std::string& path) {
    return load_npz(path).has("isect_normals");
}

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

static DualSubgridPipelineResult run_dual_for_input(const TestInput& input, size_t res) {
    if (input.is_npz())
        return run_dual_subgrid_pipeline_npz(input.npz_path);
    auto handler = input.make_handler();
    return run_dual_subgrid_pipeline(*handler, res);
}

static std::vector<TestInput> all_dual_inputs() {
    std::vector<TestInput> inputs;
    for (const char* sdf : {"Sphere", "Girl", "Cables", "Teapot", "UprightPiano", "Cube", "GrandPiano"})
        inputs.push_back({sdf, [=]{ return std::make_unique<SDFQueryHandler>(sdf, 1e-2f); }, ""});
    for (const auto& path : find_files_with_ext(TEST_DATA_DIR, {".obj", ".ply", ".stl", ".off"}))
        inputs.push_back({fs::path(path).stem().string(),
                          [=]{ return make_mesh_handler(path); }, ""});
    for (const auto& path : find_files_with_ext(TEST_NPZ_DIR, {".npz"}))
        if (npz_has_normals(path)) // dual needs intersection normals
            inputs.push_back({"npz:" + fs::path(path).stem().string(), nullptr, path});
    return inputs;
}

// Expand inputs into (input, res) cases: sweep resolutions for SDF/mesh inputs,
// run each npz once (it carries its own resolution).
static std::vector<DualCase> all_dual_cases() {
    std::vector<DualCase> cases;
    for (const auto& in : all_dual_inputs()) {
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

// ── tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("Dual Subgrid pipeline: primal mesh is manifold and edge-manifold",
          "[dual_subgrid_pipeline]") {
    auto cases = all_dual_cases();
    if (cases.size() <= 1)
        WARN("Only SDF inputs available — set TEST_DATA_DIR for mesh coverage");

    auto tc = GENERATE_COPY(from_range(cases));
    const auto& input = tc.input;
    auto res = tc.res;
    INFO("input: " << input.name << (input.is_npz() ? "  (npz, own resolution)"
                                                    : "  res: " + std::to_string(res)));

    auto result = run_dual_for_input(input, res);

    if (result.soup.faces.empty())
        SKIP("No polygons produced");

    auto [mesh, geo] = construct_primal_mesh_from_face_per_edge_data(
        result.soup.faces, result.soup.vertices, result.soup.faces_per_edge
    );
    MeshOwner owner{mesh, geo};

    REQUIRE(mesh != nullptr);
    CHECK(mesh->isManifold());
    CHECK(check_edge_manifoldness(*mesh));
    INFO("non_even_tets: " << result.non_even_tets);
    if (result.non_even_tets == 0) {
        mesh->greedilyOrientFaces();
        CHECK(mesh->isOriented());
    }
}

TEST_CASE("Dual Subgrid pipeline: dual mesh can be constructed",
          "[dual_subgrid_pipeline]") {
    auto tc = GENERATE_COPY(from_range(all_dual_cases()));
    const auto& input = tc.input;
    auto res = tc.res;
    INFO("input: " << input.name << (input.is_npz() ? "  (npz, own resolution)"
                                                    : "  res: " + std::to_string(res)));

    auto result = run_dual_for_input(input, res);

    if (result.soup.faces.empty())
        SKIP("No polygons produced");

    auto [primal_mesh, primal_geo] = construct_primal_mesh_from_face_per_edge_data(
        result.soup.faces, result.soup.vertices, result.soup.faces_per_edge
    );
    MeshOwner primal_owner{primal_mesh, primal_geo};

    REQUIRE(primal_mesh->isManifold());

    auto [dual_mesh, dual_geo] = construct_dual_mesh(
        *primal_mesh, *primal_geo, result.soup.dual_positions
    );
    MeshOwner dual_owner{dual_mesh, dual_geo};

    REQUIRE(dual_mesh != nullptr);
    CHECK(dual_mesh->nFaces() > 0);
    CHECK(dual_mesh->nVertices() > 0);
}
