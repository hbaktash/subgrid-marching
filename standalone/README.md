# Standalone subgrid ports (copy-pasteable, zero dependencies)

Self-contained ports of the **per-tet** Subgrid Marching Tetrahedra local
constructions — primal and dual — with **no external dependencies**. Drop a
single file into your project and call one of two functions.

- **C++:** [`cpp/subgrid_mt.hpp`](cpp/subgrid_mt.hpp) — header-only, STL only
  (the dual QEF solve is a hand-rolled 3×3, so no Eigen).
- **Python:** [`python/subgrid_mt.py`](python/subgrid_mt.py) — standard library
  only (NumPy is **not** required).

These reproduce the maintained core in `src/subgrid_MT/` **exactly** (verified
in-process against it; the only numerical difference is ~1e-11 in the dual QEF
point, from the Eigen→hand-rolled solver swap).

## Quick use

Both entry points take one tetrahedron's data — 4 vertex positions, 4 global
vertex indices, and 6 edges' worth of intersection parameters `t ∈ [0,1]`
(edge order `(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)`) — matching
[`docs/construction_policy.md`](../docs/construction_policy.md) → "Input".
`boundary_comb_curves` is folded inside, so there's nothing else to call.

- `subgrid_primal(positions, indices, edge_isect_ts)` → primal surface
- `subgrid_dual(positions, indices, edge_isect_ts, edge_isect_normals, reg_alpha, project_duals)`
  → dual (QEF) surface; needs one surface normal per intersection

Each returns a plain struct/dataclass with the locally-merged geometry
(`vertices`, `faces`, …) plus a per-vertex combinatorial **signature**. To build
a global mesh across many tets, merge vertices whose `signature_key(...)` match
(a `None`/empty key means "never merge — allocate fresh"). See `examples/` for a
ready-to-read global-merge snippet.

### C++

```bash
cd examples
c++ -std=c++17 -I../cpp primal_demo.cpp -o primal_demo && ./primal_demo
c++ -std=c++17 -I../cpp dual_demo.cpp   -o dual_demo   && ./dual_demo
```

### Python

```bash
cd examples
python3 primal_demo.py
python3 dual_demo.py
```

## What each thing is for

| Path | Purpose | Ships with your copy? |
|------|---------|-----------------------|
| `cpp/subgrid_mt.hpp` | the C++ port | **yes — this is the product** |
| `python/subgrid_mt.py` | the Python port | **yes — this is the product** |
| `examples/` | tiny primal+dual demos incl. a cross-tet global-merge snippet | no (reference only) |
| `verify/` | harnesses that check the ports against the repo core (C++ in-process + fixture-based for Python) | no |
| `CMakeLists.txt` | self-contained project that builds **only** the verification harnesses | no |

If you just want the algorithm, copy `cpp/subgrid_mt.hpp` **or**
`python/subgrid_mt.py` and ignore everything else.

## Verifying

The ports are checked against the repository core (the ground truth), not against
each other. `standalone/` is its own CMake project (it is **not** wired into the
repo-root build); it pulls this repo in only to build the core library it checks
against, so all artifacts land in `standalone/build/`.

```bash
# Configure standalone/ as its own project (headless: no OpenGL/GLFW needed).
cmake -S standalone -B standalone/build -DSUBGRID_POLYSCOPE_VIEWER=OFF
cmake --build standalone/build --target verify_port verify_port_dual

# C++ ports, checked in-process against the core (prints PASS / mismatch counts).
./standalone/build/verify_port
./standalone/build/verify_port_dual

# Non-C++ ports: dump the core's reference fixtures, then check the port against them.
./standalone/build/verify_port      standalone/verify/fixtures_primal.txt
python3 standalone/verify/verify_primal.py standalone/verify/fixtures_primal.txt
```

The C++ harnesses enumerate all edge-int configs up to S4 symmetry plus random
configs, and compare vertices, faces, and signatures exactly (dual points to
~1e-6, from the Eigen→hand-rolled QEF swap). The fixture `.txt` files are large
and git-ignored — regenerate them whenever the core or a harness changes.
