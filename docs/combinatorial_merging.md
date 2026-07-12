# Vertex merging

Per-tet reconstruction produces a triangle *soup*: adjacent tetrahedra each emit
their own copy of the vertices they share. The three policies below all emit the
**same triangles** — they differ only in how (or whether) those shared duplicate
vertices are welded. `subgrid` exposes all three; the dual pipeline always uses
combinatorial merge.

## Policies at a glance

| Policy | Flag | Output | Speed |
|--------|------|--------|-------|
| **Combinatorial** (default) | *(none)* | Edge-manifold, and watertight + orientable when every tet is even-sum. Exact and resolution-independent — no tolerance to tune. | Slowest (adds ~16–19% construction time) |
| **Numerical** | `--mergeEPS <eps>` | Positional welding of vertices within `<eps>` (KD-tree union-find). Best-effort, prone to error. | Faster |
| **None** | `--noMerge` | Raw soup — one vertex per tet corner, with duplicates at every shared face. | Fastest |

Merging only removes duplicate *vertices*, never triangles, so all three produce
the identical face count. Combinatorial merge reduces the vertex count about 4–5×
versus the raw soup.

## Which one to use

- **Combinatorial merge is the default and the recommended choice.** It is exact
  and resolution-independent, so the result is edge-manifold — and watertight and
  orientable whenever every tet is even-sum (see the note on non-even tets in the
  [README](../README.md)). It costs a bit of construction time. The dual pipeline
  always uses this mode.
- **Numerical merge** is faster but error-prone. A single global epsilon cannot both stitch every shared-face
  duplicate (adjacent tets recompute them with small floating-point drift) *and*
  avoid collapsing genuinely-distinct nearby vertices, so on detailed inputs it
  can leave cracks (eps too small) or non-manifold pinches (eps too large). Prefer
  combinatorial merge unless you specifically need positional welding.
- **No merge** is the fastest and is handy for inspecting the raw per-tet output,
  but leaves duplicate vertices along every shared face. Could be used for positional welding as a post-process via other programs.

---

## How combinatorial merging works

Instead of comparing floating-point positions, each output vertex carries a
**combinatorial signature** (a `CombVertexKey`) that describes *where* it lives in
the grid. Two vertices are merged iff their signatures are equal, so the result is
exact and independent of resolution or intersection spacing.

### The vertex signature

A signature has the fields `(i, j, order, k, type)`:

- **`i`, `j`** — the endpoints of the tet edge the vertex sits on, normalized so
  `i < j`. They are tet-local (0–3) during per-tet construction and converted to
  global grid-node ids before the global merge.
- **`order`** — the vertex's 1-based position along edge `{i, j}`. For a regular
  vertex this is its intersection order; for a Steiner vertex it is the lower of
  the two consecutive positions it sits between.
- **`k`** — the third face vertex, used only by face-shared Steiner vertices;
  `-1` otherwise.
- **`type`** — determines the *scope* over which the vertex is allowed to merge
  (see below). An extra `edge_int` field (total intersections on `{i, j}`) is
  stored but is **not** part of the identity; it is only used to flip `order` when
  normalizing the edge direction.

### Merge scope by type

| Type | Merges across | Key |
|------|---------------|-----|
| `NORMAL` | any tet sharing edge `{i, j}` | `(i, j, order)` |
| `SCOOP_FACE_STEINER` | tets sharing face `{i, j, k}` | `(i, j, order, k)` |
| `SCOOP_INTERIOR_STEINER` | within a single tet only | `(i, j, order)` |
| `STANDALONE` | never — always a fresh vertex | — (e.g. octagon centers, spiral means, recursion centroids) |

Interior Steiner vertices only dedup *within* a tet because their position depends
on the tet's interior geometry, which differs per tet; standalone vertices never
merge at all.

### Normalization and remapping

Signatures are always stored with `i < j`. When a raw vertex comes in with
`i > j`, its endpoints and `order` are flipped (regular vertices flip as
`order = edge_int − order + 1`; Steiner vertices as `order = edge_int − order`).
As per-tet soups are assembled into the global mesh, indices are remapped from
tet-local to global node ids and orders are converted from edge-local to
global-within-tet, keeping signatures consistent so shared vertices collide
exactly. Getting this flip/remap logic right is the crux of the exact merge.

The dual pipeline uses the same machinery but only ever produces `NORMAL`
vertices, so it needs none of the Steiner/standalone bookkeeping.
