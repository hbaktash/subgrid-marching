# Performance: threading and query reuse

Three options control how the tet loop runs. Only threading is opt-in; the other
two are on by default and exist mainly to be turned off. All compose, and none
changes the surface the per-tet construction produces.

| Option | CLI | Python | Default |
|--------|-----|--------|---------|
| Worker threads | `-j, --threads` | `num_threads=` | `1` (serial); `0` = all cores |
| Edge query reuse | `--noQueryCache` to disable | `query_cache=` | **on** |
| Canonical query direction | `--noCanonicalQueries` to disable | `canonical_queries=` | **on** |

```sh
./build/subgrid -i mesh.obj -r 128 -j 0 -o out.obj      # caching is already on
```

```python
smt.primal_from_mesh_file("mesh.obj", 128, num_threads=0)
```

Measured on a 10-core M-series laptop, primal extraction unless noted:

| run | serial, uncached | `-j 0` (cache on by default) | |
|---|---|---|---|
| `-s Cables -r 64` (SDF) | 15.9s | **1.39s** | 11.5× |
| `-i spot.obj -r 128` (mesh) | 7.21s | **2.57s** | 2.8× |
| `dualSubgrid -i spot.obj -r 128` | 6.68s | **3.61s** | 1.85× |
| `-i spot.obj -r 48 --cgal` (exact) | 2.93s | **0.42s** | 7.0× |

---

## Why there is anything to save

Reconstruction is per-tetrahedron, but the *inputs* are not. In the 5-tet
decomposition of each grid cube, an edge is shared by 4 tets if it is
axis-aligned and 6 if it is a face diagonal. The plain loop asks for all six
edges of every tet, so over an `n³` grid it issues `5n³ × 6 = 30n³` edge queries
against only

```
3n(n+1)² + 3n²(n+1)  ≈  6n³
```

distinct edges — **about five queries per edge**. Meanwhile the per-tet work
itself is embarrassingly parallel; only the final merge is not.

## Edge query reuse (on by default)

The cache reaches exactly one query per unique edge (asserted in
`tests/test_pipeline_parity.cpp`). A tet in cube layer `iz` only touches grid
nodes in planes `iz` and `iz+1`, so the cache holds a two-plane window and slides
it. The window's contents follow from index arithmetic alone, with no tet
iteration and no queries, which has two consequences: memory is O(n²) rather than
O(n³), and the fill phase can pre-size its output so threads write disjoint slots
with no locking. Because plane `iz+1` becomes the next window's plane `iz`,
nothing is ever queried twice.

Only edges that actually have crossings are stored — inside a loaded plane a
lookup miss unambiguously means "queried, none", because the plane was enumerated
exhaustively.

| input | `--noQueryCache` | default |
|---|---|---|
| mesh `spot.obj` r=128 | query 3.02s | **1.71s** |
| SDF `Cables` r=64 | query 13.1s | **2.81s** |
| `--cgal` `spot.obj` r=48 | total 2.93s | **1.13s** |

**`slab` needs the implicit grid, and that is the only place a cache helps.**
For `--npz` input it is ignored: `PrecomputedQueryHandler` is
*already* a hash lookup on the canonical edge.

## Canonical query direction (on by default)

The cache asks for each edge once, in the canonical `min(i,j) → max(i,j)`
direction, and flips the result for tets that traverse it the other way. The
uncached path instead asks once per incident tet, in that tet's own direction —
and **no handler is exactly direction-symmetric**; either due to nature of the query in sphere marching or numerical errors in ray triangle intersection (except with CGAL EPECK).

Canonicalizing makes the construction produce the same mesh regardless of how each tet vertices are locally or globally indexed, and regardless of which end of the edge the ray starts from. Tested canonicalization against the uncached path on a bunch of examples:

- The output for meshes does not generally change with/without canonicalization, because the mesh query is already symmetric. The only exception is when the mesh goes through a grid node, which is a degenerate case; the flags `--seed` (perturbing the input) and `--cgal` (Exact queries) are exactly made to take care of these type of degeneracies.
- For SDFs, canonicalization is important: given how we march the SDF and jump across the roots, the output can change depending on which end of the edge the ray starts from.  
  

## `-j`: threads

The tet loop runs in chunks — one cube layer per chunk for the grid. Within a
chunk the query and per-tet construction run across threads; the chunk is then
merged into the global soup **sequentially, in flat tet order**.

Threading uses `std::thread` only (`include/common/parallel_for.h`) — no OpenMP,
no TBB. Chunks are handed out dynamically from an atomic cursor rather than split
into fixed per-thread ranges, because per-tet cost varies by orders of magnitude
(an empty tet versus one with several boundary curves) and the non-empty tets
cluster near the surface.

### What does not scale

The merge is sequential, so it bounds the speedup: it is roughly 10-20% of the tet
loop, capping it near 5-6× however many cores you have. The dual pipeline gains
least (1.85×) because its post-pipeline assembly —
`construct_primal_mesh_from_face_per_edge_data` and `construct_dual_mesh` — is
sequential and untouched by this work; at high resolution that, not the tet loop,
is its bottleneck.

### Timing fields under threading

`isect_time` and `construction_time` split as their names suggest only in serial
runs. With threads the query and construction phases overlap, so summing
per-thread times would exceed elapsed: `isect_time` then holds the wall-clock
time of the parallel phases and `construction_time` the sequential merge. They
still add up to the tet loop's elapsed time.

## Thread safety

The per-tet code in `subgrid_MT/` has no mutable shared state and is untouched by
this work. Of the query handlers:

- **FCPW** (`MeshQueryHandler`) — `fcpw::Scene::intersect` is `const`; safe.
- **SDF** — the built-in registry is a read-only `static const` map; safe.
- **Precomputed** — reads its edge map, writes only out-parameters; safe.
- **CGAL** (`CGALQueryHandler`) — safe, but *only* because the build defines
  `CGAL_HAS_THREADS`. EPECK's lazy kernel reference-counts its handles, and those
  counters are non-atomic without it. CGAL's `config.h` derives the macro from
  `BOOST_HAS_THREADS` / `_OPENMP` — neither of which is set here, since the tet
  loop uses `std::thread` — so `CMakeLists.txt` defines it explicitly under
  `SUBGRID_WITH_CGAL=ON`. Without that define the handler would race silently.
  This is the one place where getting it wrong is a correctness bug rather than a
  slowdown, so the parity and thread tests include a CGAL input specifically to
  put it under ThreadSanitizer.

  `--cgal` is also where threading pays best (7×): exact predicates make each
  query expensive enough that the sequential merge barely registers.

The handler's diagnostic counters (`total_queries`, `edge_query_counts`) are
plain non-atomic shared state; the threaded path turns them off via
`collect_query_stats` and skips `update_global_query_count_map`.
