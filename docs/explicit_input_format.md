# Explicit Tet-Mesh + Precomputed Intersections Input Format

This document specifies the `.npz` input format for providing an explicit
tetrahedral mesh together with precomputed edge–surface intersections. This
allows external code (e.g. a GPU pipeline) to compute intersection data and
feed it directly into the subgrid/dual-subgrid reconstruction pipelines, bypassing the
built-in mesh/SDF query handlers.

## File format: NumPy `.npz`

A single `.npz` archive with the following named arrays:

| Array name       | Shape              | Dtype   | Description |
|------------------|--------------------|---------|-------------|
| `vertices`       | `(V, 3)`           | float64 | Tet-mesh vertex positions |
| `tets`           | `(T, 4)`           | int32   | Tet connectivity (indices into `vertices`) |
| `edges`          | `(E, 2)`           | int32   | Edges with ≥1 intersection (indices into `vertices`) |
| `isect_offsets`  | `(E+1,)`           | int32   | CSR-style offsets into packed intersection arrays |
| `isect_ts`       | `(N_total,)`       | float64 | Packed t-values for all intersections |
| `isect_normals`  | `(N_total, 3)`     | float64 | *(optional)* Surface normals at each intersection |

Where `N_total = isect_offsets[E]` is the total number of intersection points
across all edges.

## Conventions

### Edge orientation and t-value parameterization

Edges in the `edges` array MUST be stored with the smaller vertex index first:

    edges[k] = (i, j)  where  i < j

The t-values for edge `(i, j)` parameterize the segment from vertex `i` to
vertex `j`:

    P(t) = (1 - t) * vertices[i]  +  t * vertices[j]

so `t = 0` corresponds to vertex `i` and `t = 1` to vertex `j`.

### Sorted t-values

Within each edge's segment `isect_ts[offsets[k] : offsets[k+1]]`, the t-values
MUST be sorted in ascending order:

    t_1 ≤ t_2 ≤ ... ≤ t_m

### Reverse-edge lookup

When the reconstruction pipeline processes a tet whose local edge runs from
global vertex `j` to global vertex `i` (with `j > i`), the lookup finds the
stored entry under `(i, j)` and reverses the intersection data:

- t-values become: `1 - t_m, 1 - t_{m-1}, ..., 1 - t_1`
- normals (if present) are reordered to match: `n_m, n_{m-1}, ..., n_1`

Normals are surface normals and do NOT get negated — only their order reverses
to stay associated with the correct (flipped) t-values.

### Accessing intersections for an edge

For edge index `k` (0-based), the intersections are:

    ts     = isect_ts[isect_offsets[k] : isect_offsets[k+1]]
    norms  = isect_normals[isect_offsets[k] : isect_offsets[k+1]]  // if present

### Vertex indexing and mesh stitching

The tet mesh SHOULD use shared vertex indices for vertices at coincident
positions (i.e., a conforming tet mesh, not a "tet soup"). This enables
combinatorial merging of the per-tet triangle soups: two tets sharing an edge
will look up the same `(i, j)` key, get identical t-values, and produce
boundary vertices that match exactly without numerical tolerance.

If the input is a tet soup (duplicated vertices at shared faces), the
reconstruction still works but the output mesh can only be stitched via
numerical epsilon-merging (`--mergeEPS`). The user is responsible for this
choice.

## Python writer example

```python
import numpy as np

def save_explicit_input(path, vertices, tets, edge_pairs, edge_ts, edge_normals=None):
    """
    Parameters
    ----------
    vertices : (V, 3) float64
    tets : (T, 4) int32
    edge_pairs : list of (i, j) with i < j, length E
    edge_ts : list of 1D arrays, one sorted array of t-values per edge
    edge_normals : list of (m_k, 3) arrays, or None
    """
    edges = np.array(edge_pairs, dtype=np.int32)
    offsets = np.zeros(len(edge_ts) + 1, dtype=np.int32)
    for k, ts in enumerate(edge_ts):
        offsets[k + 1] = offsets[k] + len(ts)

    packed_ts = np.concatenate(edge_ts) if edge_ts else np.empty(0, dtype=np.float64)

    arrays = {
        'vertices': np.asarray(vertices, dtype=np.float64),
        'tets': np.asarray(tets, dtype=np.int32),
        'edges': edges,
        'isect_offsets': offsets,
        'isect_ts': packed_ts,
    }
    if edge_normals is not None:
        packed_normals = np.concatenate(edge_normals) if edge_normals else np.empty((0, 3), dtype=np.float64)
        arrays['isect_normals'] = packed_normals

    np.savez(path, **arrays)
```

## CLI usage

```bash
# Primal subgrid pipeline (combinatorial merge is the default)
./subgrid --npz input.npz -o output.obj

# Dual subgrid pipeline
./dualSubgrid --npz input.npz -o output.obj
```