// bindings.cpp — pybind11 bindings for Subgrid Marching Tetrahedra.
//
// Exposes two layers of the maintained core (src/subgrid_MT/, src/*_pipeline.cpp):
//
//   1. Per-tet constructions — subgrid_primal / subgrid_dual / subgrid_greedy.
//      These fold boundary_comb_curves() into the construction call exactly the way
//      standalone/cpp/subgrid_mt.hpp does, so the Python API matches the ports in
//      standalone/ while calling the real core. Input/output follow
//      docs/construction_policy.md.
//
//   2. Full pipelines — primal_from_mesh / _from_mesh_file / _from_sdf / _from_npz
//      and their dual counterparts, mirroring the `subgrid` and `dualSubgrid` CLIs
//      including input preprocessing and output mesh assembly.
//
// Array policy: every returned buffer is a NumPy array that either owns moved C++
// memory through a capsule (vertices, dual points, signatures — no element copies)
// or is built once in a single pass (the CSR face arrays, which have no contiguous
// C++ counterpart to borrow). Inputs are read straight out of NumPy buffers and
// marshalled into the core's std::array<std::vector<...>,6> containers; that copy is
// forced by the core signatures and is bounded by the size of one tet's data.
//
// NOTE: TriangleSoup::COMB_MERGE is a process-wide static in the core, so these
// functions are not safe to call concurrently from multiple threads.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "common/triangle_soup.h"
#include "common/utils.h"
#include "nc/normalCoord_containers.h"
#include "subgrid_MT/boundary_curve.h"
#include "subgrid_MT/dual/dual_construction.h"
#include "subgrid_MT/primal_reconstruction.h"

#include "assembly/dual_mesh_assembly.h"
#include "assembly/mesh_processing.h"
#include "dual_subgrid_pipeline.h"
#include "query/cgal_queries.h"
#include "query/input_query_handler.h"
#include "subgrid_pipeline.h"

#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

// ===========================================================================
// Layout guarantees that make the zero-copy returns legal
// ===========================================================================

static_assert(std::is_standard_layout<Vector3>::value, "Vector3 must be standard layout");
static_assert(sizeof(Vector3) == 3 * sizeof(double), "Vector3 must be 3 packed doubles");
static_assert(offsetof(Vector3, x) == 0 * sizeof(double), "Vector3 layout");
static_assert(offsetof(Vector3, y) == 1 * sizeof(double), "Vector3 layout");
static_assert(offsetof(Vector3, z) == 2 * sizeof(double), "Vector3 layout");

static_assert(std::is_standard_layout<CombVertexKey>::value, "CombVertexKey must be standard layout");
static_assert(sizeof(int) == sizeof(int32_t), "CombVertexKey columns assume 32-bit int");
static_assert(sizeof(CombVertexKey) == 6 * sizeof(int32_t), "CombVertexKey must be 6 packed ints");
static_assert(offsetof(CombVertexKey, i) == 0 * sizeof(int32_t), "signature column 0 is i");
static_assert(offsetof(CombVertexKey, j) == 1 * sizeof(int32_t), "signature column 1 is j");
static_assert(offsetof(CombVertexKey, order) == 2 * sizeof(int32_t), "signature column 2 is order");
static_assert(offsetof(CombVertexKey, edge_int) == 3 * sizeof(int32_t), "signature column 3 is edge_int");
static_assert(offsetof(CombVertexKey, k) == 4 * sizeof(int32_t), "signature column 4 is k");
static_assert(offsetof(CombVertexKey, type) == 5 * sizeof(int32_t), "signature column 5 is type");

// ===========================================================================
// Zero-copy array construction: hand a moved std::vector to NumPy, keeping it
// alive with a capsule that deletes it when the last view goes away.
// ===========================================================================

template <typename T>
py::capsule make_owner(std::vector<T>* held) {
    return py::capsule(held, [](void* p) { delete static_cast<std::vector<T>*>(p); });
}

// (N, 3) float64 over a vector<Vector3>.
py::array vec3_array(std::vector<Vector3>&& v) {
    const auto n = static_cast<py::ssize_t>(v.size());
    if (n == 0) return py::array_t<double>(std::vector<py::ssize_t>{0, 3});
    auto* held = new std::vector<Vector3>(std::move(v));
    return py::array_t<double>(
        {n, py::ssize_t(3)},
        {py::ssize_t(sizeof(Vector3)), py::ssize_t(sizeof(double))},
        reinterpret_cast<const double*>(held->data()),
        make_owner(held));
}

// (N, 6) int32 over a vector<CombVertexKey>: columns [i, j, order, edge_int, k, type].
py::array signature_array(std::vector<CombVertexKey>&& v) {
    const auto n = static_cast<py::ssize_t>(v.size());
    if (n == 0) return py::array_t<int32_t>(std::vector<py::ssize_t>{0, 6});
    auto* held = new std::vector<CombVertexKey>(std::move(v));
    return py::array_t<int32_t>(
        {n, py::ssize_t(6)},
        {py::ssize_t(sizeof(CombVertexKey)), py::ssize_t(sizeof(int32_t))},
        reinterpret_cast<const int32_t*>(held->data()),
        make_owner(held));
}

// (N,) int64 over a vector<int64_t>.
py::array int64_array(std::vector<int64_t>&& v) {
    const auto n = static_cast<py::ssize_t>(v.size());
    if (n == 0) return py::array_t<int64_t>(std::vector<py::ssize_t>{0});
    auto* held = new std::vector<int64_t>(std::move(v));
    return py::array_t<int64_t>({n}, {py::ssize_t(sizeof(int64_t))}, held->data(), make_owner(held));
}

// (N, 3) int64 over a flat vector<int64_t> of triplets.
py::array int64_triplet_array(std::vector<int64_t>&& v) {
    const auto n = static_cast<py::ssize_t>(v.size() / 3);
    if (n == 0) return py::array_t<int64_t>(std::vector<py::ssize_t>{0, 3});
    auto* held = new std::vector<int64_t>(std::move(v));
    return py::array_t<int64_t>(
        {n, py::ssize_t(3)},
        {py::ssize_t(3 * sizeof(int64_t)), py::ssize_t(sizeof(int64_t))},
        held->data(),
        make_owner(held));
}

// ===========================================================================
// Result types
// ===========================================================================

// Faces are general n-gons, so they are returned in CSR form: face f spans
// face_vertices[face_offsets[f] : face_offsets[f + 1]].
struct SoupArrays {
    py::array vertices;       // (N, 3) float64
    py::array face_offsets;   // (F + 1,) int64
    py::array face_vertices;  // (M,) int64
    py::array signatures;     // (N, 6) int32, or (0, 6) when not applicable
};

struct PrimalTetResult : SoupArrays {
    bool non_even = false;
};

struct DualTetResult : SoupArrays {
    py::array dual_positions;       // (F, 3) float64
    py::array face_edge_tet_faces;  // (M, 3) int64
    bool non_even = false;
    bool non_normal = false;
};

struct PipelineResult : SoupArrays {
    int non_zero_tets = 0;
    int non_normal_tets = 0;
    int non_even_tets = 0;
    size_t total_tets = 0;
    double isect_time = 0.0;
    double construction_time = 0.0;
    double assembly_time = 0.0;
};

struct DualPipelineResult : PipelineResult {
    py::array dual_positions;       // (F, 3) float64, empty when assembled
    py::array face_edge_tet_faces;  // (M, 3) int64, empty when assembled
    bool assembled = true;
};

py::ssize_t num_faces(const SoupArrays& s) {
    const auto n = s.face_offsets.size();
    return n == 0 ? 0 : static_cast<py::ssize_t>(n) - 1;
}

// A list of NumPy views into face_vertices — slicing, so no data is copied.
py::list faces_as_views(const SoupArrays& s) {
    py::list out;
    auto offsets = py::cast<py::array_t<int64_t>>(s.face_offsets);
    auto r = offsets.unchecked<1>();
    for (py::ssize_t f = 0; f + 1 < r.shape(0); ++f)
        out.append(s.face_vertices[py::slice(static_cast<py::ssize_t>(r(f)),
                                             static_cast<py::ssize_t>(r(f + 1)), 1)]);
    return out;
}

SoupArrays soup_to_arrays(TriangleSoup&& soup) {
    SoupArrays out;

    const bool sigs_parallel = soup.signatures.size() == soup.vertices.size();

    size_t total_face_verts = 0;
    for (const auto& f : soup.faces) total_face_verts += f.size();

    std::vector<int64_t> offsets;
    offsets.reserve(soup.faces.size() + 1);
    offsets.push_back(0);
    std::vector<int64_t> flat;
    flat.reserve(total_face_verts);
    for (const auto& f : soup.faces) {
        for (size_t idx : f) flat.push_back(static_cast<int64_t>(idx));
        offsets.push_back(static_cast<int64_t>(flat.size()));
    }

    out.vertices = vec3_array(std::move(soup.vertices));
    out.face_offsets = int64_array(std::move(offsets));
    out.face_vertices = int64_array(std::move(flat));
    out.signatures = sigs_parallel ? signature_array(std::move(soup.signatures))
                                   : py::array(py::array_t<int32_t>(std::vector<py::ssize_t>{0, 6}));
    return out;
}

// faces_per_edge is parallel to faces, one tet-face triplet per polygon edge.
py::array faces_per_edge_array(const TriangleSoup& soup) {
    if (soup.faces_per_edge.empty()) return py::array_t<int64_t>(std::vector<py::ssize_t>{0, 3});
    if (soup.faces_per_edge.size() != soup.faces.size())
        throw std::runtime_error("faces_per_edge is not parallel to faces");

    std::vector<int64_t> flat;
    size_t total = 0;
    for (const auto& fpe : soup.faces_per_edge) total += fpe.size();
    flat.reserve(3 * total);
    for (size_t f = 0; f < soup.faces_per_edge.size(); ++f) {
        if (soup.faces_per_edge[f].size() != soup.faces[f].size())
            throw std::runtime_error("faces_per_edge[f] must have one entry per edge of face f");
        for (const auto& t : soup.faces_per_edge[f]) {
            flat.push_back(static_cast<int64_t>(t[0]));
            flat.push_back(static_cast<int64_t>(t[1]));
            flat.push_back(static_cast<int64_t>(t[2]));
        }
    }
    return int64_triplet_array(std::move(flat));
}

// ===========================================================================
// Input parsing
// ===========================================================================

std::string shape_str(const py::array& a) {
    std::ostringstream os;
    os << "(";
    for (py::ssize_t d = 0; d < a.ndim(); ++d) os << (d ? ", " : "") << a.shape(d);
    os << ")";
    return os.str();
}

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IntArray = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;

std::array<Vector3, 4> parse_tet_positions(const py::object& obj) {
    auto arr = DoubleArray::ensure(obj);
    if (!arr) throw py::type_error("tet_positions: expected a float array-like of shape (4, 3)");
    if (arr.ndim() != 2 || arr.shape(0) != 4 || arr.shape(1) != 3)
        throw py::value_error("tet_positions: expected shape (4, 3), got " + shape_str(arr));
    const double* p = arr.data();
    std::array<Vector3, 4> out;
    for (int v = 0; v < 4; ++v) out[v] = Vector3{p[3 * v + 0], p[3 * v + 1], p[3 * v + 2]};
    return out;
}

std::array<size_t, 4> parse_tet_indices(const py::object& obj) {
    auto arr = IntArray::ensure(obj);
    if (!arr) throw py::type_error("tet_global_indices: expected an integer array-like of shape (4,)");
    if (arr.ndim() != 1 || arr.shape(0) != 4)
        throw py::value_error("tet_global_indices: expected shape (4,), got " + shape_str(arr));
    const int64_t* p = arr.data();
    std::array<size_t, 4> out;
    for (int v = 0; v < 4; ++v) {
        if (p[v] < 0) throw py::value_error("tet_global_indices: indices must be non-negative");
        out[v] = static_cast<size_t>(p[v]);
    }
    return out;
}

// Per edge, the sorted intersection parameters. Edge order is
// (0,1), (0,2), (0,3), (1,2), (1,3), (2,3); t runs from the lower to the higher
// local vertex, so it must be non-decreasing within each edge.
std::array<std::vector<double>, 6> parse_edge_ts(const py::object& obj) {
    if (py::len(obj) != 6)
        throw py::value_error("edge_isect_ts: expected 6 entries (one per tet edge), got " +
                              std::to_string(py::len(obj)));
    std::array<std::vector<double>, 6> out;
    for (int e = 0; e < 6; ++e) {
        auto arr = DoubleArray::ensure(obj[py::int_(e)]);
        if (!arr)
            throw py::type_error("edge_isect_ts[" + std::to_string(e) +
                                 "]: expected a 1-D float array-like of t-values");
        if (arr.ndim() != 1)
            throw py::value_error("edge_isect_ts[" + std::to_string(e) + "]: expected a 1-D array, got " +
                                  shape_str(arr));
        const double* p = arr.data();
        out[e].assign(p, p + arr.size());
        for (size_t k = 0; k < out[e].size(); ++k) {
            if (!std::isfinite(out[e][k]))
                throw py::value_error("edge_isect_ts[" + std::to_string(e) + "]: t-values must be finite");
            if (k > 0 && out[e][k] < out[e][k - 1])
                throw py::value_error("edge_isect_ts[" + std::to_string(e) +
                                      "]: t-values must be sorted ascending along the edge "
                                      "(from the lower to the higher local vertex index)");
        }
    }
    return out;
}

std::array<std::vector<Vector3>, 6> parse_edge_normals(const py::object& obj,
                                                       const std::array<std::vector<double>, 6>& ts) {
    std::array<std::vector<Vector3>, 6> out;
    if (obj.is_none()) return out;
    if (py::len(obj) != 6)
        throw py::value_error("edge_isect_normals: expected 6 entries (one per tet edge), got " +
                              std::to_string(py::len(obj)));
    for (int e = 0; e < 6; ++e) {
        auto arr = DoubleArray::ensure(obj[py::int_(e)]);
        if (!arr)
            throw py::type_error("edge_isect_normals[" + std::to_string(e) +
                                 "]: expected a float array-like of shape (n, 3)");
        const auto n = static_cast<size_t>(ts[e].size());
        if (n == 0 && arr.size() == 0) continue;
        if (arr.ndim() != 2 || arr.shape(1) != 3)
            throw py::value_error("edge_isect_normals[" + std::to_string(e) +
                                  "]: expected shape (n, 3), got " + shape_str(arr));
        if (static_cast<size_t>(arr.shape(0)) != n)
            throw py::value_error("edge_isect_normals[" + std::to_string(e) + "]: has " +
                                  std::to_string(arr.shape(0)) + " normals but edge_isect_ts[" +
                                  std::to_string(e) + "] has " + std::to_string(n) + " intersections");
        const double* p = arr.data();
        out[e].resize(n);
        for (size_t k = 0; k < n; ++k) out[e][k] = Vector3{p[3 * k + 0], p[3 * k + 1], p[3 * k + 2]};
    }
    return out;
}

std::vector<Vector3> parse_positions(const py::object& obj, const char* name) {
    auto arr = DoubleArray::ensure(obj);
    if (!arr) throw py::type_error(std::string(name) + ": expected a float array-like of shape (V, 3)");
    if (arr.ndim() != 2 || arr.shape(1) != 3)
        throw py::value_error(std::string(name) + ": expected shape (V, 3), got " + shape_str(arr));
    const double* p = arr.data();
    std::vector<Vector3> out(static_cast<size_t>(arr.shape(0)));
    for (size_t v = 0; v < out.size(); ++v) out[v] = Vector3{p[3 * v + 0], p[3 * v + 1], p[3 * v + 2]};
    return out;
}

// Accepts an (F, K) integer array of same-sized polygons, or any sequence of
// per-face index sequences (ragged is fine).
std::vector<std::vector<size_t>> parse_polygons(const py::object& obj, size_t n_vertices, const char* name) {
    std::vector<std::vector<size_t>> out;

    auto check = [&](int64_t idx) {
        if (idx < 0 || static_cast<size_t>(idx) >= n_vertices)
            throw py::value_error(std::string(name) + ": vertex index " + std::to_string(idx) +
                                  " out of range for " + std::to_string(n_vertices) + " vertices");
        return static_cast<size_t>(idx);
    };

    if (py::isinstance<py::array>(obj)) {
        auto arr = IntArray::ensure(obj);
        if (!arr) throw py::type_error(std::string(name) + ": expected an integer array of shape (F, K)");
        if (arr.ndim() != 2 || arr.shape(1) < 3)
            throw py::value_error(std::string(name) + ": expected shape (F, K) with K >= 3, got " +
                                  shape_str(arr));
        const int64_t* p = arr.data();
        const auto k = static_cast<size_t>(arr.shape(1));
        out.resize(static_cast<size_t>(arr.shape(0)));
        for (size_t f = 0; f < out.size(); ++f) {
            out[f].resize(k);
            for (size_t c = 0; c < k; ++c) out[f][c] = check(p[f * k + c]);
        }
        return out;
    }

    for (const auto& face : obj) {
        std::vector<size_t> poly;
        for (const auto& idx : face) poly.push_back(check(py::cast<int64_t>(idx)));
        if (poly.size() < 3)
            throw py::value_error(std::string(name) + ": every face needs at least 3 vertices");
        out.push_back(std::move(poly));
    }
    return out;
}

// ===========================================================================
// Per-tet constructions (docs/construction_policy.md)
// ===========================================================================

std::array<int, 6> counts_from(const std::array<std::vector<double>, 6>& ts) {
    std::array<int, 6> counts{};
    for (int e = 0; e < 6; ++e) counts[e] = static_cast<int>(ts[e].size());
    return counts;
}

PrimalTetResult py_subgrid_primal(const py::object& positions, const py::object& indices,
                                  const py::object& ts, bool scoop_mid_vertices, double scoop_bulge) {
    const auto tet_positions = parse_tet_positions(positions);
    const auto tet_indices = parse_tet_indices(indices);
    const auto edge_ts = parse_edge_ts(ts);

    TriangleSoup::COMB_MERGE = true;
    EdgeOccupations non_normal_edge_occupations;
    auto curves = boundary_comb_curves(tet_indices, counts_from(edge_ts), non_normal_edge_occupations,
                                       /*populate_normal_curves=*/false);
    const auto& open_curves = std::get<0>(curves);

    TriangleSoup soup = subgrid_surface(tet_positions, tet_indices, edge_ts, open_curves, std::get<1>(curves),
                                        non_normal_edge_occupations, scoop_mid_vertices, scoop_bulge,
                                        /*verbose=*/false);

    PrimalTetResult out;
    static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(soup));
    out.non_even = !open_curves.empty();
    return out;
}

PrimalTetResult py_subgrid_greedy(const py::object& positions, const py::object& indices,
                                  const py::object& ts, bool scoop_mid_vertices, double scoop_bulge) {
    const auto tet_positions = parse_tet_positions(positions);
    const auto tet_indices = parse_tet_indices(indices);
    const auto edge_ts = parse_edge_ts(ts);

    TriangleSoup::COMB_MERGE = true;
    EdgeOccupations occupations;
    auto curves = boundary_comb_curves(tet_indices, counts_from(edge_ts), occupations,
                                       /*populate_normal_curves=*/true);
    const auto& open_curves = std::get<0>(curves);

    TriangleSoup soup = subgrid_greedy(tet_positions, tet_indices, edge_ts, std::get<1>(curves),
                                       std::get<2>(curves), scoop_mid_vertices, scoop_bulge);

    PrimalTetResult out;
    static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(soup));
    out.non_even = !open_curves.empty();
    return out;
}

DualTetResult py_subgrid_dual(const py::object& positions, const py::object& indices, const py::object& ts,
                              const py::object& normals, double reg_alpha, bool project_duals,
                              bool use_normals) {
    const auto tet_positions = parse_tet_positions(positions);
    const auto tet_indices = parse_tet_indices(indices);
    const auto edge_ts = parse_edge_ts(ts);
    if (use_normals && normals.is_none())
        throw py::value_error(
            "subgrid_dual: edge_isect_normals is required for the QEF solve. Pass normals, or "
            "use_normals=False to place each dual point at its boundary-polygon centroid.");
    const auto edge_normals = parse_edge_normals(normals, edge_ts);

    TriangleSoup::COMB_MERGE = true;
    EdgeOccupations occupations;
    auto curves = boundary_comb_curves(tet_indices, counts_from(edge_ts), occupations,
                                       /*populate_normal_curves=*/true);
    const auto& open_curves = std::get<0>(curves);
    const auto& scoop_curves = std::get<1>(curves);

    TriangleSoup soup = dual_subgrid_surface(tet_positions, tet_indices, edge_ts, edge_normals, open_curves,
                                             scoop_curves, std::get<2>(curves), reg_alpha, project_duals,
                                             use_normals);

    DualTetResult out;
    out.face_edge_tet_faces = faces_per_edge_array(soup);
    out.dual_positions = vec3_array(std::move(soup.dual_positions));
    static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(soup));
    out.non_even = !open_curves.empty();
    out.non_normal = !(open_curves.empty() && scoop_curves.empty());
    return out;
}

// ===========================================================================
// Pipeline plumbing
// ===========================================================================

enum class Merge { COMBINATORIAL, NUMERICAL, NONE };

Merge parse_merge(const std::string& mode) {
    if (mode == "combinatorial") return Merge::COMBINATORIAL;
    if (mode == "numerical") return Merge::NUMERICAL;
    if (mode == "none") return Merge::NONE;
    throw py::value_error("merge: expected 'combinatorial', 'numerical', or 'none', got '" + mode + "'");
}

std::unique_ptr<InputQueryHandler> make_mesh_handler(std::vector<Vector3> positions,
                                                     std::vector<std::vector<size_t>> polygons,
                                                     bool preprocess, unsigned int seed, bool cgal) {
    if (preprocess) {
        SimplePolygonMesh simple(polygons, positions);
        auto pre = preprocess_input_mesh(simple, /*normalize_scale=*/0.98, seed);
        positions = std::move(pre.positions);
        polygons = std::move(pre.polygons);
    }
    if (cgal) {
#ifdef HAVE_CGAL
        return std::make_unique<CGALQueryHandler>(positions, polygons);
#else
        throw std::runtime_error(
            "cgal=True requires building with -DSUBGRID_WITH_CGAL=ON (see the README's "
            "'Robust queries (CGAL)' section).");
#endif
    }
    return std::make_unique<MeshQueryHandler>(positions, polygons);
}

std::unique_ptr<InputQueryHandler> make_sdf_handler(const std::string& name, double step_size) {
    if (!SDFQueryHandler::is_valid_sdf_name(name))
        throw py::value_error("unknown SDF '" + name + "'; call available_sdfs() for the list");
    return std::make_unique<SDFQueryHandler>(name, static_cast<float>(step_size));
}

std::pair<std::vector<Vector3>, std::vector<std::vector<size_t>>> read_mesh_file(const std::string& path) {
    SimplePolygonMesh simple;
    simple.readMeshFromFile(path, "");
    return {simple.vertexCoordinates, simple.polygons};
}

std::vector<Vector3> mesh_positions(SurfaceMesh& mesh, VertexPositionGeometry& geo) {
    std::vector<Vector3> out;
    out.reserve(mesh.nVertices());
    for (Vertex v : mesh.vertices()) out.push_back(geo.inputVertexPositions[v]);
    return out;
}

// Edge intersection reuse (see query/edge_isect_cache.h). "slab" needs the
// implicit grid and is ignored for explicit tet meshes, whose intersections are
// already precomputed per edge.
QueryCache parse_query_cache(const std::string& name) {
    if (name == "none") return QueryCache::NONE;
    if (name == "slab") return QueryCache::SLAB;
    throw std::invalid_argument("query_cache must be 'none' or 'slab'; got '" + name + "'");
}

SubgridPipelineOpts make_primal_opts(bool mod2, bool greedy, bool robust, bool progress, bool verbose,
                                     double scoop_bulge, bool scoop_mid_vertices,
                                     const std::string& query_cache, int num_threads,
                                     bool canonical_queries) {
    SubgridPipelineOpts opts;
    opts.query_cache = parse_query_cache(query_cache);
    opts.num_threads = num_threads;
    opts.canonical_queries = canonical_queries;
    opts.mod2 = mod2;
    opts.greedy = greedy;
    opts.use_robust = robust;
    opts.show_progress = progress;
    opts.verbose = verbose;
    opts.scoop_bulge = scoop_bulge;
    opts.scoop_mid_vertices = scoop_mid_vertices;
    return opts;
}

DualSubgridPipelineOpts make_dual_opts(bool mod2, bool robust, bool progress, double reg_alpha,
                                       bool project_duals, bool use_normals,
                                       const std::string& query_cache, int num_threads,
                                       bool canonical_queries) {
    DualSubgridPipelineOpts opts;
    opts.query_cache = parse_query_cache(query_cache);
    opts.num_threads = num_threads;
    opts.canonical_queries = canonical_queries;
    opts.mod2 = mod2;
    opts.use_robust = robust;
    opts.show_progress = progress;
    opts.reg_alpha = reg_alpha;
    opts.project_duals = project_duals;
    opts.no_normal = !use_normals;
    return opts;
}

// Mirrors the postprocessing in apps/main_subgrid.cpp: merge policy, then (for the
// combinatorial default) a greedy orientation pass over the assembled mesh.
PipelineResult primal_finish(SubgridPipelineResult&& r, Merge merge, double merge_eps, bool orient) {
    PipelineResult out;
    out.non_zero_tets = r.non_zero_tets;
    out.non_normal_tets = r.non_normal_tets;
    out.non_even_tets = r.non_even_tets;
    out.total_tets = r.total_tets;
    out.isect_time = r.isect_time;
    out.construction_time = r.construction_time;

    TriangleSoup assembled;
    const auto t_start = std::chrono::high_resolution_clock::now();
    {
        py::gil_scoped_release nogil;
        if (merge == Merge::NONE) {
            assembled = std::move(r.soup);
        } else if (merge == Merge::NUMERICAL) {
            std::unique_ptr<SurfaceMesh> mesh;
            std::unique_ptr<VertexPositionGeometry> geo;
            auto pair = mergeIdenticalVertices(merge_eps, r.soup.faces, r.soup.vertices);
            mesh.reset(pair.first);
            geo.reset(pair.second);
            mesh->compress();
            assembled.vertices = mesh_positions(*mesh, *geo);
            assembled.faces = mesh->getFaceVertexList();
        } else {
            auto mesh_and_geo = makeSurfaceMeshAndGeometry(r.soup.faces, r.soup.vertices);
            auto& mesh = std::get<0>(mesh_and_geo);
            auto& geo = std::get<1>(mesh_and_geo);
            if (orient) mesh->greedilyOrientFaces();
            mesh->compress();
            assembled.vertices = mesh_positions(*mesh, *geo);
            assembled.faces = mesh->getFaceVertexList();
            // Signatures stay parallel as long as no vertex was dropped.
            if (assembled.vertices.size() == r.soup.signatures.size())
                assembled.signatures = std::move(r.soup.signatures);
        }
    }
    out.assembly_time =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();
    static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(assembled));
    return out;
}

// Mirrors apps/main_dual_subgrid.cpp: build the primal mesh from the face-per-edge
// tags, then the dual polygon mesh from the QEF points.
DualPipelineResult dual_finish(DualSubgridPipelineResult&& r, bool assemble) {
    DualPipelineResult out;
    out.non_zero_tets = r.non_zero_tets;
    out.non_normal_tets = r.non_normal_tets;
    out.non_even_tets = r.non_even_tets;
    out.total_tets = r.total_tets;
    out.isect_time = r.isect_time;
    out.construction_time = r.construction_time;
    out.assembled = assemble;

    if (!assemble) {
        out.face_edge_tet_faces = faces_per_edge_array(r.soup);
        out.dual_positions = vec3_array(std::move(r.soup.dual_positions));
        static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(r.soup));
        return out;
    }

    TriangleSoup assembled;
    const auto t_start = std::chrono::high_resolution_clock::now();
    {
        py::gil_scoped_release nogil;
        auto primal = construct_primal_mesh_from_face_per_edge_data(r.soup.faces, r.soup.vertices,
                                                                    r.soup.faces_per_edge);
        std::unique_ptr<SurfaceMesh> primal_mesh(primal.first);
        std::unique_ptr<VertexPositionGeometry> primal_geo(primal.second);
        if (!primal_mesh->isManifold())
            throw std::runtime_error(
                "the assembled primal mesh is not manifold, so the dual cannot be built. "
                "Pass assemble=False to get the raw per-tet boundary loops and dual points.");

        auto dual = construct_dual_mesh(*primal_mesh, *primal_geo, r.soup.dual_positions);
        std::unique_ptr<SurfaceMesh> dual_mesh(dual.first);
        std::unique_ptr<VertexPositionGeometry> dual_geo(dual.second);
        dual_mesh->compress();
        assembled.vertices = mesh_positions(*dual_mesh, *dual_geo);
        assembled.faces = dual_mesh->getFaceVertexList();
    }
    out.assembly_time =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();

    out.dual_positions = py::array_t<double>(std::vector<py::ssize_t>{0, 3});
    out.face_edge_tet_faces = py::array_t<int64_t>(std::vector<py::ssize_t>{0, 3});
    static_cast<SoupArrays&>(out) = soup_to_arrays(std::move(assembled));
    return out;
}

PipelineResult run_primal(std::unique_ptr<InputQueryHandler> handler, size_t resolution,
                          const SubgridPipelineOpts& opts, Merge merge, double merge_eps, bool orient) {
    TriangleSoup::COMB_MERGE = (merge == Merge::COMBINATORIAL);
    SubgridPipelineResult r;
    {
        py::gil_scoped_release nogil;
        r = run_subgrid_pipeline(*handler, resolution, opts);
    }
    return primal_finish(std::move(r), merge, merge_eps, orient);
}

DualPipelineResult run_dual(std::unique_ptr<InputQueryHandler> handler, size_t resolution,
                            const DualSubgridPipelineOpts& opts, bool assemble) {
    DualSubgridPipelineResult r;
    {
        py::gil_scoped_release nogil;
        r = run_dual_subgrid_pipeline(*handler, resolution, opts);
    }
    return dual_finish(std::move(r), assemble);
}

}  // namespace

// ===========================================================================
// Module
// ===========================================================================

PYBIND11_MODULE(_core, m) {
    m.doc() =
        "Subgrid Marching Tetrahedra — intersection-free isosurface extraction.\n"
        "\n"
        "Two layers are exposed:\n"
        "  * per-tet constructions (subgrid_primal / subgrid_dual / subgrid_greedy),\n"
        "    matching the reference ports in standalone/;\n"
        "  * full pipelines (primal_from_* / dual_from_*), matching the `subgrid` and\n"
        "    `dualSubgrid` command line tools.\n"
        "\n"
        "Faces are general n-gons, so they come back in CSR form (face_offsets plus\n"
        "face_vertices) with a `.faces` list of NumPy views over the same buffer.\n"
        "\n"
        "Not thread-safe: the core keeps a process-wide merge flag, so calls must not\n"
        "overlap across threads.";

    py::enum_<CombVertexSigType>(m, "CombVertexSigType", py::arithmetic(),
                                 "Kind of vertex a combinatorial signature refers to. STANDALONE "
                                 "vertices must never be merged across tets — give each a fresh index.")
        .value("NORMAL", CombVertexSigType::NORMAL)
        .value("SCOOP_FACE_STEINER", CombVertexSigType::SCOOP_FACE_STEINER)
        .value("SCOOP_INTERIOR_STEINER", CombVertexSigType::SCOOP_INTERIOR_STEINER)
        .value("STANDALONE", CombVertexSigType::STANDALONE);

    // Column indices into the (N, 6) `signatures` array.
    m.attr("SIG_I") = 0;
    m.attr("SIG_J") = 1;
    m.attr("SIG_ORDER") = 2;
    m.attr("SIG_EDGE_INT") = 3;  // not part of the vertex identity
    m.attr("SIG_K") = 4;
    m.attr("SIG_TYPE") = 5;

    {
        std::vector<int64_t> pairs;
        pairs.reserve(12);
        for (const auto& p : ALL_TET_PAIRS) {
            pairs.push_back(static_cast<int64_t>(p.first));
            pairs.push_back(static_cast<int64_t>(p.second));
        }
        auto arr = py::array_t<int64_t>(std::vector<py::ssize_t>{6, 2});
        std::memcpy(arr.mutable_data(), pairs.data(), pairs.size() * sizeof(int64_t));
        m.attr("EDGE_VERTEX_PAIRS") = arr;
    }

    m.def("edge_pair_to_index", &edge_pair_to_index, "i"_a, "j"_a,
          "Index of tet edge {i, j} in the canonical order (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).");

    m.def("available_sdfs", &SDFQueryHandler::available_sdf_names,
          "Names of the built-in signed distance functions accepted by primal_from_sdf / dual_from_sdf.");

    // ---- result types -----------------------------------------------------

    py::class_<SoupArrays>(m, "Soup",
                           "Polygon soup in CSR form. Face f is\n"
                           "face_vertices[face_offsets[f] : face_offsets[f + 1]].")
        .def_readonly("vertices", &SoupArrays::vertices, "(N, 3) float64 vertex positions.")
        .def_readonly("face_offsets", &SoupArrays::face_offsets, "(F + 1,) int64 CSR offsets into face_vertices.")
        .def_readonly("face_vertices", &SoupArrays::face_vertices, "(M,) int64 concatenated face vertex indices.")
        .def_readonly("signatures", &SoupArrays::signatures,
                      "(N, 6) int32 combinatorial signatures parallel to `vertices`: columns\n"
                      "[i, j, order, edge_int, k, type]. Two vertices from different tets are the\n"
                      "same point iff columns [i, j, order, k, type] match and type is not\n"
                      "STANDALONE. Empty when merging is off or the vertices were re-indexed.")
        .def_property_readonly("faces", &faces_as_views,
                               "Faces as a list of NumPy views into `face_vertices` (no copy).\n"
                               "Materializes one view per face, so prefer the CSR arrays for large outputs.")
        .def_property_readonly("num_faces", &num_faces, "Number of faces.")
        .def("__repr__", [](const SoupArrays& s) {
            return "<Soup " + std::to_string(s.vertices.shape(0)) + " vertices, " +
                   std::to_string(num_faces(s)) + " faces>";
        });

    py::class_<PrimalTetResult, SoupArrays>(m, "PrimalResult", "Result of a single-tet primal construction.")
        .def_readonly("non_even", &PrimalTetResult::non_even,
                      "True when this tet is not even-sum: the boundary curves do not close up, so\n"
                      "the surface will be open here. See the README's 'non-even (open) tets' note.")
        .def("__repr__", [](const PrimalTetResult& r) {
            return "<PrimalResult " + std::to_string(r.vertices.shape(0)) + " vertices, " +
                   std::to_string(num_faces(r)) + " faces, non_even=" + (r.non_even ? "True" : "False") + ">";
        });

    py::class_<DualTetResult, SoupArrays>(m, "DualResult", "Result of a single-tet dual construction.")
        .def_readonly("dual_positions", &DualTetResult::dual_positions,
                      "(F, 3) float64, one QEF dual point per boundary loop, parallel to the faces.")
        .def_readonly("face_edge_tet_faces", &DualTetResult::face_edge_tet_faces,
                      "(M, 3) int64, parallel to `face_vertices`: the tet face {i, j, k} each polygon\n"
                      "edge lies on. Edge e of face f is the segment from face_vertices[e] to the next\n"
                      "vertex of that face, wrapping at the end. Used to build the global twin map.")
        .def_readonly("non_even", &DualTetResult::non_even, "True when the tet is not even-sum (open curves).")
        .def_readonly("non_normal", &DualTetResult::non_normal, "True when the tet has open or scoop curves.")
        .def("__repr__", [](const DualTetResult& r) {
            return "<DualResult " + std::to_string(r.vertices.shape(0)) + " vertices, " +
                   std::to_string(num_faces(r)) + " loops, non_even=" + (r.non_even ? "True" : "False") + ">";
        });

    py::class_<PipelineResult, SoupArrays>(m, "PipelineResult", "Result of a full primal extraction.")
        .def_readonly("non_zero_tets", &PipelineResult::non_zero_tets, "Tets the surface passed through.")
        .def_readonly("non_normal_tets", &PipelineResult::non_normal_tets,
                      "Active tets with open or scoop curves.")
        .def_readonly("non_even_tets", &PipelineResult::non_even_tets,
                      "Active tets that failed the even-sum condition. Zero on a watertight input;\n"
                      "a positive count means the output may have holes or pinch vertices.")
        .def_readonly("total_tets", &PipelineResult::total_tets, "Tets in the grid.")
        .def_readonly("isect_time", &PipelineResult::isect_time, "Seconds spent in intersection queries.")
        .def_readonly("construction_time", &PipelineResult::construction_time,
                      "Seconds spent in the per-tet construction.")
        .def_readonly("assembly_time", &PipelineResult::assembly_time,
                      "Seconds spent merging and orienting the output mesh.")
        .def("__repr__", [](const PipelineResult& r) {
            return "<PipelineResult " + std::to_string(r.vertices.shape(0)) + " vertices, " +
                   std::to_string(num_faces(r)) + " faces, non_even_tets=" +
                   std::to_string(r.non_even_tets) + ">";
        });

    py::class_<DualPipelineResult, PipelineResult>(m, "DualPipelineResult",
                                                   "Result of a full dual extraction.")
        .def_readonly("assembled", &DualPipelineResult::assembled,
                      "True when `vertices` / `faces` hold the assembled dual mesh; False when they\n"
                      "hold the raw per-tet boundary loops.")
        .def_readonly("dual_positions", &DualPipelineResult::dual_positions,
                      "(F, 3) float64 dual point per boundary loop. Only filled when assemble=False.")
        .def_readonly("face_edge_tet_faces", &DualPipelineResult::face_edge_tet_faces,
                      "(M, 3) int64 tet face per polygon edge. Only filled when assemble=False.")
        .def("__repr__", [](const DualPipelineResult& r) {
            return "<DualPipelineResult " + std::to_string(r.vertices.shape(0)) + " vertices, " +
                   std::to_string(num_faces(r)) + (r.assembled ? " dual faces" : " boundary loops") +
                   ", non_even_tets=" + std::to_string(r.non_even_tets) + ">";
        });

    // ---- per-tet constructions -------------------------------------------

    m.def("subgrid_primal", &py_subgrid_primal, "tet_positions"_a, "tet_global_indices"_a, "edge_isect_ts"_a,
          py::kw_only(), "scoop_mid_vertices"_a = true, "scoop_bulge"_a = 0.001,
          R"doc(Primal (even-sum) reconstruction inside a single tetrahedron.

Args:
    tet_positions: (4, 3) float64 — the tet's vertex positions.
    tet_global_indices: (4,) int — the vertices' ids in the global tet mesh. Adjacent
        tets must pass consistent ids: they drive the combinatorial signatures that
        make shared boundary points merge exactly.
    edge_isect_ts: 6 sequences of intersection parameters, one per edge in the order
        (0,1), (0,2), (0,3), (1,2), (1,3), (2,3). For edge (i, j) with i < j the value
        t denotes (1 - t) * v_i + t * v_j, and values must be sorted ascending.
    scoop_mid_vertices: insert simplicial-embedding Steiner vertices (default on).
    scoop_bulge: push-in epsilon for those Steiner vertices. 1e-3 matches the
        pipeline default documented in the README; standalone/ uses 1e-4.

Returns:
    PrimalResult with `vertices`, CSR faces, per-vertex `signatures`, and `non_even`.

Input normals are not used by the primal construction.)doc");

    m.def("subgrid_greedy", &py_subgrid_greedy, "tet_positions"_a, "tet_global_indices"_a, "edge_isect_ts"_a,
          py::kw_only(), "scoop_mid_vertices"_a = true, "scoop_bulge"_a = 0.001,
          R"doc(Alternative primal construction: fan / spiral triangulation of each boundary curve.

Same arguments and result as subgrid_primal. Faster and simpler, but unlike
subgrid_primal the output is not guaranteed to be free of self-intersections.)doc");

    m.def("subgrid_dual", &py_subgrid_dual, "tet_positions"_a, "tet_global_indices"_a, "edge_isect_ts"_a,
          "edge_isect_normals"_a = py::none(), py::kw_only(), "reg_alpha"_a = 0.1,
          "project_duals"_a = false, "use_normals"_a = true,
          R"doc(Dual (QEF / dual-contouring) reconstruction inside a single tetrahedron.

Args:
    tet_positions, tet_global_indices, edge_isect_ts: as in subgrid_primal.
    edge_isect_normals: 6 arrays of shape (n_e, 3) — a surface normal per intersection,
        parallel to edge_isect_ts. Required unless use_normals=False.
    reg_alpha: QEF regularization; larger values pull dual points toward the centroid
        of their boundary polygon.
    project_duals: clip each dual point back inside the tet.
    use_normals: when False, normals are ignored and each dual point is placed at its
        boundary-polygon centroid — no QEF solve, no normals needed.

Returns:
    DualResult: the boundary loops (`vertices`, CSR faces), one `dual_positions` point
    per loop, and `face_edge_tet_faces` tagging the tet face of every polygon edge.)doc");

    // ---- primal pipelines -------------------------------------------------

    m.def(
        "primal_from_mesh",
        [](const py::object& vertices, const py::object& faces, size_t resolution, bool preprocess,
           unsigned int seed, bool cgal, bool mod2, bool greedy, double scoop_bulge, bool scoop_mid_vertices,
           const std::string& merge, double merge_eps, bool orient, bool robust, bool progress, bool verbose,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            auto V = parse_positions(vertices, "vertices");
            auto F = parse_polygons(faces, V.size(), "faces");
            const Merge merge_mode = parse_merge(merge);
            const auto opts = make_primal_opts(mod2, greedy, robust, progress, verbose, scoop_bulge,
                                               scoop_mid_vertices, query_cache, num_threads, canonical_queries);
            std::unique_ptr<InputQueryHandler> handler;
            {
                py::gil_scoped_release nogil;
                handler = make_mesh_handler(std::move(V), std::move(F), preprocess, seed, cgal);
            }
            return run_primal(std::move(handler), resolution, opts, merge_mode, merge_eps, orient);
        },
        "vertices"_a, "faces"_a, "resolution"_a = 64, py::kw_only(), "preprocess"_a = true, "seed"_a = 1u,
        "cgal"_a = false, "mod2"_a = false, "greedy"_a = false, "scoop_bulge"_a = 0.001,
        "scoop_mid_vertices"_a = true, "merge"_a = "combinatorial", "merge_eps"_a = 1e-8, "orient"_a = true,
        "robust"_a = false, "progress"_a = false, "verbose"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        R"doc(Extract a primal isosurface from an in-memory triangle mesh.

Args:
    vertices: (V, 3) float64 mesh vertex positions.
    faces: (F, K) integer array of same-sized polygons, or a sequence of index
        sequences. Polygons are triangulated during preprocessing.
    resolution: grid resolution r; builds an r^3 cube grid (5 r^3 tetrahedra).
    preprocess: weld, triangulate, normalize into 98% of [-1, 1]^3, and apply the
        fixed-seed rigid offset that decorrelates the mesh from the grid. Leave this
        on unless your mesh is already inside the grid box.
    seed: seed for that rigid offset; change it to shake off degenerate hits.
    cgal: use the exact CGAL (EPECK) query handler. Requires -DSUBGRID_WITH_CGAL=ON.
    mod2: reduce each edge's intersection count mod 2 before construction.
    greedy: use the greedy fan construction instead of the even-sum one.
    merge: "combinatorial" (exact, default), "numerical" (positional, uses merge_eps),
        or "none" (raw per-tet soup).
    merge_eps: welding radius for merge="numerical". No single value is right for
        every input: too small leaves cracks, too large creates pinches. The default
        only welds bitwise-coincident vertices, which still leaves the per-tet scoop
        Steiner vertices apart. Prefer the combinatorial default.
    orient: greedily orient the faces of the merged mesh.
    robust: robust ray-intersection queries (roughly 2x slower, rarely needed).
    progress: print a progress bar to the process stdout.

Returns:
    PipelineResult. Check `non_even_tets`: it should be 0 for a watertight input.)doc");

    m.def(
        "primal_from_mesh_file",
        [](const std::string& path, size_t resolution, bool preprocess, unsigned int seed, bool cgal,
           bool mod2, bool greedy, double scoop_bulge, bool scoop_mid_vertices, const std::string& merge,
           double merge_eps, bool orient, bool robust, bool progress, bool verbose,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const Merge merge_mode = parse_merge(merge);
            const auto opts = make_primal_opts(mod2, greedy, robust, progress, verbose, scoop_bulge,
                                               scoop_mid_vertices, query_cache, num_threads, canonical_queries);
            std::unique_ptr<InputQueryHandler> handler;
            {
                py::gil_scoped_release nogil;
                auto mesh = read_mesh_file(path);
                handler = make_mesh_handler(std::move(mesh.first), std::move(mesh.second), preprocess, seed,
                                            cgal);
            }
            return run_primal(std::move(handler), resolution, opts, merge_mode, merge_eps, orient);
        },
        "path"_a, "resolution"_a = 64, py::kw_only(), "preprocess"_a = true, "seed"_a = 1u, "cgal"_a = false,
        "mod2"_a = false, "greedy"_a = false, "scoop_bulge"_a = 0.001, "scoop_mid_vertices"_a = true,
        "merge"_a = "combinatorial", "merge_eps"_a = 1e-8, "orient"_a = true, "robust"_a = false,
        "progress"_a = false, "verbose"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a primal isosurface from a mesh file (OBJ / PLY / OFF).\n"
        "Same options as primal_from_mesh; this is the exact equivalent of\n"
        "`subgrid -i <path> -r <resolution>`.");

    m.def(
        "primal_from_sdf",
        [](const std::string& name, size_t resolution, double step_size, bool mod2, bool greedy,
           double scoop_bulge, bool scoop_mid_vertices, const std::string& merge, double merge_eps,
           bool orient, bool robust, bool progress, bool verbose,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const Merge merge_mode = parse_merge(merge);
            const auto opts = make_primal_opts(mod2, greedy, robust, progress, verbose, scoop_bulge,
                                               scoop_mid_vertices, query_cache, num_threads, canonical_queries);
            auto handler = make_sdf_handler(name, step_size);
            return run_primal(std::move(handler), resolution, opts, merge_mode, merge_eps, orient);
        },
        "name"_a, "resolution"_a = 64, py::kw_only(), "step_size"_a = 1e-3, "mod2"_a = false,
        "greedy"_a = false, "scoop_bulge"_a = 0.001, "scoop_mid_vertices"_a = true,
        "merge"_a = "combinatorial", "merge_eps"_a = 1e-8, "orient"_a = true, "robust"_a = false,
        "progress"_a = false, "verbose"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a primal isosurface from a built-in signed distance function.\n"
        "`name` must be one of available_sdfs(); `step_size` is the marching step used\n"
        "to bracket sign changes along each tet edge.");

    m.def(
        "primal_from_npz",
        [](const std::string& path, bool mod2, bool greedy, double scoop_bulge, bool scoop_mid_vertices,
           const std::string& merge, double merge_eps, bool orient, bool progress, bool verbose,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const Merge merge_mode = parse_merge(merge);
            const auto opts = make_primal_opts(mod2, greedy, /*robust=*/false, progress, verbose, scoop_bulge,
                                               scoop_mid_vertices, query_cache, num_threads, canonical_queries);
            TriangleSoup::COMB_MERGE = (merge_mode == Merge::COMBINATORIAL);
            SubgridPipelineResult r;
            {
                py::gil_scoped_release nogil;
                r = run_subgrid_pipeline_npz(path, opts);
            }
            return primal_finish(std::move(r), merge_mode, merge_eps, orient);
        },
        "path"_a, py::kw_only(), "mod2"_a = false, "greedy"_a = false, "scoop_bulge"_a = 0.001,
        "scoop_mid_vertices"_a = true, "merge"_a = "combinatorial", "merge_eps"_a = 1e-8, "orient"_a = true,
        "progress"_a = false, "verbose"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a primal isosurface from an explicit tet mesh with precomputed edge\n"
        "intersections (.npz). The tet mesh comes from the file, so there is no\n"
        "resolution argument. See docs/explicit_input_format.md for the layout.");

    // ---- dual pipelines ---------------------------------------------------

    m.def(
        "dual_from_mesh",
        [](const py::object& vertices, const py::object& faces, size_t resolution, bool preprocess,
           unsigned int seed, bool cgal, bool mod2, double reg_alpha, bool project_duals, bool use_normals,
           bool assemble, bool robust, bool progress,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            auto V = parse_positions(vertices, "vertices");
            auto F = parse_polygons(faces, V.size(), "faces");
            const auto opts = make_dual_opts(mod2, robust, progress, reg_alpha, project_duals, use_normals, query_cache, num_threads, canonical_queries);
            std::unique_ptr<InputQueryHandler> handler;
            {
                py::gil_scoped_release nogil;
                handler = make_mesh_handler(std::move(V), std::move(F), preprocess, seed, cgal);
            }
            return run_dual(std::move(handler), resolution, opts, assemble);
        },
        "vertices"_a, "faces"_a, "resolution"_a = 64, py::kw_only(), "preprocess"_a = true, "seed"_a = 1u,
        "cgal"_a = false, "mod2"_a = false, "reg_alpha"_a = 0.1, "project_duals"_a = false,
        "use_normals"_a = true, "assemble"_a = true, "robust"_a = false, "progress"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        R"doc(Extract a dual (QEF) surface from an in-memory triangle mesh.

Shares the mesh arguments and preprocessing of primal_from_mesh. Dual-specific:
    reg_alpha: QEF regularization weight.
    project_duals: clip each dual point back inside its grid cell.
    use_normals: when False, dual points go to boundary-polygon centroids (no QEF).
    assemble: when True (default) the result holds the assembled dual polygon mesh,
        exactly what `dualSubgrid` writes. When False it holds the raw per-tet
        boundary loops plus `dual_positions` and `face_edge_tet_faces`, so you can
        run your own assembly.)doc");

    m.def(
        "dual_from_mesh_file",
        [](const std::string& path, size_t resolution, bool preprocess, unsigned int seed, bool cgal,
           bool mod2, double reg_alpha, bool project_duals, bool use_normals, bool assemble, bool robust,
           bool progress,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const auto opts = make_dual_opts(mod2, robust, progress, reg_alpha, project_duals, use_normals, query_cache, num_threads, canonical_queries);
            std::unique_ptr<InputQueryHandler> handler;
            {
                py::gil_scoped_release nogil;
                auto mesh = read_mesh_file(path);
                handler = make_mesh_handler(std::move(mesh.first), std::move(mesh.second), preprocess, seed,
                                            cgal);
            }
            return run_dual(std::move(handler), resolution, opts, assemble);
        },
        "path"_a, "resolution"_a = 64, py::kw_only(), "preprocess"_a = true, "seed"_a = 1u, "cgal"_a = false,
        "mod2"_a = false, "reg_alpha"_a = 0.1, "project_duals"_a = false, "use_normals"_a = true,
        "assemble"_a = true, "robust"_a = false, "progress"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a dual (QEF) surface from a mesh file (OBJ / PLY / OFF).\n"
        "Equivalent to `dualSubgrid -i <path> -r <resolution>`.");

    m.def(
        "dual_from_sdf",
        [](const std::string& name, size_t resolution, double step_size, bool mod2, double reg_alpha,
           bool project_duals, bool use_normals, bool assemble, bool robust, bool progress,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const auto opts = make_dual_opts(mod2, robust, progress, reg_alpha, project_duals, use_normals, query_cache, num_threads, canonical_queries);
            auto handler = make_sdf_handler(name, step_size);
            return run_dual(std::move(handler), resolution, opts, assemble);
        },
        "name"_a, "resolution"_a = 64, py::kw_only(), "step_size"_a = 1e-3, "mod2"_a = false,
        "reg_alpha"_a = 0.1, "project_duals"_a = false, "use_normals"_a = true, "assemble"_a = true,
        "robust"_a = false, "progress"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a dual (QEF) surface from a built-in signed distance function.\n"
        "`name` must be one of available_sdfs().");

    m.def(
        "dual_from_npz",
        [](const std::string& path, bool mod2, double reg_alpha, bool project_duals, bool use_normals,
           bool assemble, bool progress,
           const std::string& query_cache, int num_threads, bool canonical_queries) {
            const auto opts = make_dual_opts(mod2, /*robust=*/false, progress, reg_alpha, project_duals,
                                             use_normals, query_cache, num_threads, canonical_queries);
            DualSubgridPipelineResult r;
            {
                py::gil_scoped_release nogil;
                try {
                    r = run_dual_subgrid_pipeline_npz(path, opts);
                } catch (const std::exception& e) {
                    // The core points at the CLI flag; say it in Python terms instead.
                    std::string msg = e.what();
                    const std::string cli = "Pass --noNormal to";
                    const auto at = msg.find(cli);
                    if (at != std::string::npos) msg.replace(at, cli.size(), "Pass use_normals=False to");
                    throw std::runtime_error(msg);
                }
            }
            return dual_finish(std::move(r), assemble);
        },
        "path"_a, py::kw_only(), "mod2"_a = false, "reg_alpha"_a = 0.1, "project_duals"_a = false,
        "use_normals"_a = true, "assemble"_a = true, "progress"_a = false,
        "query_cache"_a = "none", "num_threads"_a = 1, "canonical_queries"_a = true,
        "Extract a dual (QEF) surface from an explicit tet mesh with precomputed edge\n"
        "intersections (.npz). The archive needs an `isect_normals` array unless\n"
        "use_normals=False. See docs/explicit_input_format.md.");
}
