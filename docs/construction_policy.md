# Construction Policy: Inputs and Outputs

This document specifies the input the Subgrid Marching Tetrahedra reconstruction
accepts and the output it produces, for both the primal and dual pipelines.

## Input

Reconstruction is driven **per tetrahedron**. The core inputs for a single tet
are:

- **`tet_positions`** — `array<Vector3, 4>`: the 3D coordinates of the four tet
  vertices.
- **`tet_global_indices`** — `array<size_t, 4>`: the vertices' ids in the global
  tet grid/mesh. Used for exact combinatorial merging and for
  orientation-independent tie-breaking, so that adjacent tets agree on shared
  boundary points.
- **`edge_isect_ts`** — `array<vector<double>, 6>`: per edge, the list of
  surface-intersection parameters `t` along that edge.
- **`edge_isect_normals`** — `array<vector<Vector3>, 6>` *(optional)*: per edge, a
  surface normal at each intersection point, in the same order/sizes as
  `edge_isect_ts`. Required only by the dual pipeline; the primal pipeline
  ignores it.

### Edge ordering and parameterization

The six edges are always ordered `(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)`
(`ALL_TET_PAIRS`). For edge `(i, j)` with `i < j`, a parameter `t` denotes the
point `P(t) = (1 − t)·v_i + t·v_j`, so `t = 0` is `v_i` and `t = 1` is `v_j`.
Within each edge, t-values are stored **sorted ascending in the canonical
`i < j` direction**. This convention is relied upon throughout.

### Input sources

Three interchangeable front-ends produce the per-tet edge data:

- **Triangle mesh** (OBJ/OFF/PLY) — `MeshQueryHandler`. Intersections and normals
  are found by casting a ray along each tet edge through an FCPW BVH. This needs
  no inside/outside test and no consistent surface orientation.
- **Closed-form SDF** — `SDFQueryHandler`. Marches along each edge with a minimum
  step and brackets sign changes.
- **Explicit tet mesh + precomputed intersections** (`.npz`) —
  `PrecomputedQueryHandler`. External code (e.g. a GPU pipeline) supplies
  vertices, tets, edges, and intersections directly. See
  [explicit_input_format.md](explicit_input_format.md).

## Output

All pipelines accumulate results into a `TriangleSoup`
([triangle_soup.h](../include/common/triangle_soup.h)):

- **`vertices`** — `vector<Vector3>`.
- **`faces`** — `vector<vector<size_t>>`: polygons indexing into `vertices`. Faces
  may be general n-gons, not only triangles.
- **`signatures`** — `vector<CombVertexKey>` parallel to `vertices`; populated only
  when combinatorial merging is enabled.
- **`faces_per_edge`**, **`dual_positions`** — populated only by the dual pipeline
  (see below).

### Primal construction

The primal pipeline discards input normals. Two constructions are available:

- **`subgrid_surface`** — fills each traced boundary curve with an
  intersection-free spanning disk (the normal / even-sum construction).
- **`subgrid_greedy`** — an alternative that fans each boundary curve directly
  (triangle or centroid spiral); selected with `--greedy`.

The soup is emitted in one of three merge modes, chosen at the application level:

- **Raw soup (no merge)** — vertices may be duplicated across adjacent tets.
- **Numerical merge** (`--mergeEPS`) — merges coincident vertices within an
  epsilon (KD-tree). Error-free merging yields a manifold, orientable mesh
  whenever every active tet is even-sum (reported as `non_even_tets == 0`).
- **Combinatorial merge** (default, exact) — deduplicates by combinatorial
  signature `CombVertexKey (i, j, order, k, type)` rather than by position: no
  epsilon and no duplicate vertices. Each output vertex carries the signature
  identifying which global grid edge it lies on and its order along that edge
  (interior/Steiner vertices are marked accordingly). See
  [combinatorial_merging.md](combinatorial_merging.md).

### Dual construction

The dual pipeline keeps the boundary **loops** (rather than filled disks) and
additionally records:

- **`faces_per_edge`** — for each polygon edge, the tet face `{i, j, k}` it lies
  on. These tags build a global twin map that makes the primal mesh manifold and
  lets the dual be a clean per-vertex polygon.
- **`dual_positions`** — one dual point per polygon, placed by a QEF
  ("dual contouring") solve over the polygon's intersection points and normals,
  which recovers sharp features. Regularized by `reg_alpha`; optionally clipped
  back inside the tet with `project_duals`.

## Entry points

Library ([subgrid_pipeline.h](../include/subgrid_pipeline.h),
[dual_subgrid_pipeline.h](../include/dual_subgrid_pipeline.h)):

- `run_subgrid_pipeline(handler, resolution)` /
  `run_subgrid_pipeline_npz(npz_path)` → `SubgridPipelineResult { soup, non_even_tets }`.
- `run_dual_subgrid_pipeline(handler, resolution, reg_alpha, project_duals)` /
  `run_dual_subgrid_pipeline_npz(...)` →
  `DualSubgridPipelineResult { soup, non_even_tets, non_normal_tets }`.

Command-line apps: **`subgrid`** (primal) and **`dualSubgrid`** (dual). `subgrid` takes exactly
one of `--input` (mesh), `--inputSDF`, or `--npz`, plus `-r` (grid resolution),
`--mergeEPS` / `--noMerge` (merge mode; combinatorial merge is the default),
`--mod2`, and `--greedy`.
