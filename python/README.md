# subgrid_marching — Python bindings

Python bindings for [**Subgrid Marching Tetrahedra**](https://hbaktash.github.io/projects/subgrid-marching-tetrahedra/index.html),
a method for extracting isosurfaces that are **guaranteed free of self-intersections**,
using the intersections between a surface and the edges of a tetrahedral grid.

The bindings wrap the same C++ core the `subgrid` and `dualSubgrid` command line tools
use, and reproduce their output exactly. Everything comes back as NumPy arrays.

```python
import subgrid_marching as smt

result = smt.primal_from_mesh_file("spot.obj", resolution=64)
result.vertices        # (N, 3) float64
result.faces           # list of index arrays, one per face
result.non_even_tets   # 0 on a watertight input
```

## Installing

The build compiles the C++ core, which lives partly in git submodules, so install from
a **recursive clone** rather than from a tarball:

```sh
git clone --recursive https://github.com/hbaktash/subgrid-marching.git
cd subgrid-marching
pip install .
```

You need a C++17 compiler and CMake ≥ 3.15. The first build also fetches
[fcpw](https://github.com/rohan-sawhney/fcpw) and pybind11, and takes a few minutes;
the interactive viewer, the CLIs and the C++ test suite are all skipped.

If you forgot `--recursive`, `git submodule update --init --recursive` fixes it.

For development, `pip install -e .` plus `pytest python/tests`. Optional test
dependencies come from `pip install ".[test]"`.

## Two ways in

### Per tetrahedron

The heart of the method. You supply one tet and the points where the surface crosses
its six edges; you get back the local piece of the isosurface. Use this when you have
your own tet mesh and intersection data — from a GPU pass, a solver, or your own
sampling — and want to drive the construction yourself.

```python
import numpy as np
import subgrid_marching as smt

tet_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
tet_global_indices = np.array([0, 1, 2, 3])      # ids in your global tet mesh
edge_isect_ts = [[0.3], [0.4], [0.5], [], [], []]  # per edge, sorted t values

primal = smt.subgrid_primal(tet_positions, tet_global_indices, edge_isect_ts)
dual = smt.subgrid_dual(tet_positions, tet_global_indices, edge_isect_ts, normals)
```

| Function | Returns |
|---|---|
| `subgrid_primal(positions, indices, edge_isect_ts, *, scoop_mid_vertices=True, scoop_bulge=1e-3)` | `PrimalResult` — the intersection-free spanning surface |
| `subgrid_dual(positions, indices, edge_isect_ts, edge_isect_normals=None, *, reg_alpha=0.1, project_duals=False, use_normals=True)` | `DualResult` — boundary loops plus one QEF dual point each |
| `subgrid_greedy(positions, indices, edge_isect_ts, *, scoop_mid_vertices=True, scoop_bulge=1e-3)` | `PrimalResult` — fan triangulation; faster, **not** guaranteed intersection-free |

These fold the boundary-curve tracing inside, matching the reference ports in
[`standalone/`](https://github.com/hbaktash/subgrid-marching/tree/main/standalone) —
and they are checked against those ports in `python/tests`.

### Whole pipelines

The same thing the CLIs do: build a tet grid, query the input for edge intersections,
run the per-tet construction everywhere, and assemble the output mesh.

| Primal | Dual | Input |
|---|---|---|
| `primal_from_mesh(vertices, faces, resolution=64, ...)` | `dual_from_mesh(...)` | in-memory triangle mesh |
| `primal_from_mesh_file(path, resolution=64, ...)` | `dual_from_mesh_file(...)` | OBJ / PLY / OFF file |
| `primal_from_sdf(name, resolution=64, ...)` | `dual_from_sdf(...)` | built-in SDF, see `available_sdfs()` |
| `primal_from_npz(path, ...)` | `dual_from_npz(...)` | explicit tet mesh + precomputed hits |

Shared options mirror the CLI flags: `mod2`, `greedy`, `merge` /
`merge_eps`, `preprocess`, `seed`, `cgal`, `progress`, plus `num_threads`,
`query_cache` and `canonical_queries` (below). The dual adds `reg_alpha`,
`project_duals`, `use_normals`, and `assemble`. `help(smt.primal_from_mesh)` has the
full list.

#### Going faster

```python
smt.primal_from_mesh_file("mesh.obj", 128, num_threads=0, query_cache="slab")
```

| Option | Default | Effect |
|---|---|---|
| `num_threads` | `1` | Worker threads for the tet loop; `0` uses every core |
| `query_cache` | `"none"` | `"slab"` reuses edge intersections across the tets sharing an edge (grid input only; ignored for `*_from_npz`, whose hits are already per-edge) |

Roughly 2-3x on mesh input and up to 10x+ on SDFs, where each query is more
expensive. **Neither changes the result:** output is bit-identical at any thread
count, and the cached path matches the uncached one exactly.

That exactness rests on `canonical_queries=True` (the default), which queries
each grid edge once in a fixed direction instead of once per incident tet in
that tet's own direction, so tets sharing an edge always agree about where the
surface crosses it. It costs nothing. It does shift results very slightly versus
versions predating it — mesh input keeps its topology, SDF input can move a
little more — so pass `canonical_queries=False` to reproduce the older output.

See [docs/performance.md](https://github.com/hbaktash/subgrid-marching/blob/main/docs/performance.md)
for the mechanics and the measurements.

Vertex merging defaults to the exact **combinatorial** merge, which is what makes the
output watertight without an epsilon. `merge="numerical"` welds by position within
`merge_eps` and `merge="none"` returns the raw per-tet soup; both are best-effort, and
the numerical one is [prone to cracks or pinches](https://github.com/hbaktash/subgrid-marching/blob/main/docs/combinatorial_merging.md)
depending on the epsilon you pick.

Mesh input is preprocessed exactly as the CLI does it — welded, triangulated,
normalized into 98% of `[-1, 1]³`, then rigidly offset by a fixed-seed translation to
break axis-alignment with the grid. Pass `preprocess=False` only if your mesh is
already inside the grid box.

## What you get back

Faces are general n-gons (triangles through octagons for the primal; arbitrary polygons
for the dual), so they come back in CSR form:

```python
result.vertices        # (N, 3) float64
result.face_offsets    # (F + 1,) int64 — face f spans [face_offsets[f], face_offsets[f+1])
result.face_vertices   # (M,) int64 — concatenated vertex indices
result.faces           # list of NumPy views into face_vertices (no copy)
result.num_faces       # F

np.diff(result.face_offsets)          # polygon sizes, vectorized
result.face_vertices.reshape(-1, 3)   # valid when every face is a triangle
```

Pipeline results add `non_zero_tets`, `non_normal_tets`, `non_even_tets`, `total_tets`,
and the `isect_time` / `construction_time` / `assembly_time` breakdown. **`non_even_tets`
is the number worth checking**: it should be 0 for a watertight input. A non-zero count
means some tets saw an odd number of surface crossings, so the output may have small
holes or pinch vertices there — see the [non-even tets note](https://github.com/hbaktash/subgrid-marching#a-note-on-non-even-open-tets).

Dual results carry `dual_positions` (one point per boundary loop) and
`face_edge_tet_faces` (the tet face each polygon edge lies on). By default
`dual_from_*` assembles the dual polygon mesh for you; pass `assemble=False` to get the
raw per-tet loops instead.

### Signatures, and merging across tets

Per-tet results carry a **combinatorial signature** per vertex — a `(N, 6)` int32 array
with columns `[i, j, order, edge_int, k, type]`, addressable as `smt.SIG_I`, `smt.SIG_J`,
`smt.SIG_ORDER`, `smt.SIG_EDGE_INT`, `smt.SIG_K`, `smt.SIG_TYPE`.

Two vertices from different tets are the *same point* exactly when columns
`[i, j, order, k, type]` match — `edge_int` is carried along but is not part of the
identity. This is what makes the merge exact instead of epsilon-based. The one
exception: `type == CombVertexSigType.STANDALONE` must never merge; give each such
vertex a fresh index.

`examples/global_merge.py` is a complete, 20-line implementation.

## Conventions

The six edges are always ordered `(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)`, available
as `smt.EDGE_VERTEX_PAIRS` and `smt.edge_pair_to_index(i, j)`. For edge `(i, j)` with
`i < j`, a parameter `t` denotes `(1 - t) · v_i + t · v_j`, and the t values on an edge
must be **sorted ascending in that direction**. The bindings check this and raise
`ValueError` rather than producing quiet nonsense.

`tet_global_indices` are the ids of the tet's vertices in your global tet mesh.
Neighbouring tets must pass consistent ids: they drive both the signatures and the
orientation-independent tie-breaking that makes adjacent tets agree on shared points.

Full specification: [construction policy](https://github.com/hbaktash/subgrid-marching/blob/main/docs/construction_policy.md),
[explicit input format](https://github.com/hbaktash/subgrid-marching/blob/main/docs/explicit_input_format.md),
[vertex merging](https://github.com/hbaktash/subgrid-marching/blob/main/docs/combinatorial_merging.md).

## Notes

- **Memory.** Returned arrays own the C++ buffer they were built from, handed over
  rather than copied (`result.vertices.base` is the capsule that owns it). `.faces`
  entries are slices of `face_vertices`, not copies.
- **Threading.** To parallelize one extraction, pass `num_threads`; the work happens
  in C++ with the GIL released. Do *not* instead run several pipeline calls
  concurrently from Python threads — the core keeps a process-wide merge flag, so
  overlapping calls corrupt each other. One call at a time, internally parallel.
- **Progress output.** `progress=True` writes to the process stdout from C++, so it
  shows up in a terminal but not in a notebook cell.
- **Exact queries.** `cgal=True` needs the optional CGAL build, which the published
  wheels do **not** include — compile from a clone with
  `CMAKE_ARGS="-DSUBGRID_WITH_CGAL=ON" pip install .` (Boost, GMP and MPFR must be
  installed first; see
  [Robust queries (CGAL)](https://github.com/hbaktash/subgrid-marching#robust-queries-cgal)).
  Without it, `cgal=True` raises `RuntimeError`.

## Examples

- [`examples/single_tet.py`](examples/single_tet.py) — the per-tet primal and dual calls.
- [`examples/global_merge.py`](examples/global_merge.py) — stitching per-tet soups with signatures.
- [`examples/extract_and_save.py`](examples/extract_and_save.py) — full pipeline, written out as OBJ.

## Citation

```bibtex
@article{Baktash:2026:SMT,
  author    = {Baktash, Hossein and Gillespie, Mark and Crane, Keenan},
  title     = {Subgrid Marching Tetrahedra},
  journal   = {ACM Trans. Graph.},
  volume    = {45},
  number    = {4},
  articleno = {57},
  year      = {2026},
  doi       = {10.1145/3811358}
}
```

MIT licensed, like the rest of the repository.
