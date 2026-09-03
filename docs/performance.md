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
For `--npz` input it is ignored (with a warning): `PrecomputedQueryHandler` is
*already* a hash lookup on the canonical edge, so caching it would wrap one hash
lookup in another. Measured, a general reference-counted cache made the npz query
phase about twice as slow (0.096s → 0.179s on `spot_N64`) while reaching the same
query count — which is why only the slab window survives.

## Canonical query direction (on by default)

The cache asks for each edge once, in the canonical `min(i,j) → max(i,j)`
direction, and flips the result for tets that traverse it the other way. The
uncached path instead asks once per incident tet, in that tet's own direction —
and **no handler is exactly direction-symmetric**. Measured over 400 random
segments (`tests/test_pipeline_parity.cpp`), worst `|t_fwd − (1−t_rev)|`:

| handler | t asymmetry | crossing-count disagreement |
|---|---|---|
| precomputed (`.npz`) | 0 — already canonical | none |
| CGAL / EPECK | 6.0e-15 — rounding `t` to double | none in 400 |
| mesh / FCPW | 3.3e-6 — float internals | none in 400 |
| SDF `Sphere`, step 1e-2 | 2.1e-4 | none in 400 |
| SDF `Cables`, step 1e-2 | 6.7e-3 ≈ ⅔ of a step | **10 / 400 (2.5%)** |

CGAL is the exception that proves the rule: in exact arithmetic its `t` and
`1-t` are exact rationals, and the entire residue comes from rounding each to a
double independently at the very end — a few tens of ULPs, some nine orders of
magnitude tighter than FCPW.

So without this flag, **two tets sharing a grid edge already disagree** about
where — and for SDF input, sometimes how many times — the surface crosses it.
That is pre-existing behaviour, not something the cache introduced; it is simply
why a cached run does not match an uncached one byte for byte.

Canonicalizing makes the uncached path ask the same way the cache does. It issues
**no extra queries** (measured: 6.079s vs 6.079s on `Cables` r=48, 0.422s vs
0.428s on `spot` r=64) — only a reversal on the roughly half of tet edges that run
high index to low. With it on, every configuration agrees bit for bit:

| | cached vs uncached |
|---|---|
| canonical off (`--noCanonicalQueries`) | differs (see tolerances above) |
| canonical on (default) | **identical** |

It is on by default because it is free and removes a real inconsistency, but it
does change what the extractor produces. Measured against the non-canonical path:

- **10 paper meshes at r=96:** 9 kept identical topology (max vertex shift
  5e-7 – 2.7e-6, consistent with a ~1.5e-4 error in `t` over a 0.021-long grid
  edge). The tenth, `cell_anatomy1.obj`, changed by a single vertex and face.
  Sweeping every one of its 5,391,648 unique grid edges in both directions found
  exactly **3** where the crossing *count* disagreed — and all three touch the
  same grid node, `(60,53,23)`, with the disputed crossing sitting at
  `t = 3.5e-7` or `t = 0.9999996`. In other words the input surface passes almost
  exactly through that node, and whether the grazing hit is reported depends on
  which end the ray starts from. This is the degenerate coincidence `--seed`
  exists to break: at seeds 2, 3 and 4 the same mesh has **zero** count
  mismatches. It is not a general property of mesh input.
- **All 63 built-in SDFs at r=64:** `non_even` came out **identical on every
  one** — never worse, never better, including the nine with nonzero counts.
  Vertex and face counts matched exactly on 37 of 62; among those that differed
  the median was 0.19%, and the outliers were the fractals, where the surface has
  detail far below the sampling scale and the crossing count is genuinely
  ambiguous: `Serpinski` 10.2% of vertices (but 0.06% of faces), `Mandelbulb`
  6.6% / 3.8%, `Julia` 3.7% / 0.05%.

Face counts move far less than vertex counts almost everywhere, i.e. the surface
itself is essentially unchanged and most of the delta is in how many distinct
vertices the merge produces. Pass `--noCanonicalQueries` to reproduce
pre-canonical output exactly.

## `-j`: threads

The tet loop runs in chunks — one cube layer per chunk for the grid. Within a
chunk the query and per-tet construction run across threads; the chunk is then
merged into the global soup **sequentially, in flat tet order**.

The merge has to stay sequential (it threads every vertex signature through one
shared map), and keeping it in tet order is what makes the result independent of
the thread count. Per-tet construction depends only on the tet's own positions,
global indices, and edge intersections — never on iteration history — so local
soups are order-independent and replaying them in order reproduces the serial
mesh exactly. **Output is bit-identical at any thread count**, verified across
1/2/4/8 threads for both pipelines, all cache modes, and all input families.

Chunking is a requirement, not a tuning knob: every local `TriangleSoup` carries
its own signature map, so buffering all of them before merging would run past a
gigabyte at n=128.

Threading uses `std::thread` only (`include/common/parallel_for.h`) — no OpenMP,
no TBB. Chunks are handed out dynamically from an atomic cursor rather than split
into fixed per-thread ranges, because per-tet cost varies by orders of magnitude
(an empty tet versus one with several boundary curves) and the non-empty tets
cluster near the surface.

### What does not scale

The merge is sequential, so it bounds the speedup: it is roughly 10% of the tet
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

The test suite runs clean under ThreadSanitizer.
