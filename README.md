# Subgrid Marching Tetrahedra

Implementation of [**Subgrid Marching Tetrahedra**](https://hbaktash.github.io/projects/subgrid-marching-tetrahedra/index.html), 
a method for extracting intersection-free isosurfaces using intersections with a tetrahedral grid's edges.

The two primary executables are:

| Executable | Description |
|------------|-------------|
| `subgrid` | Primal extraction — outputs a triangle soup / manifold mesh |
| `dualSubgrid` | Dual extraction — outputs a dual polygon mesh via QEF solve |

An interactive single-tetrahedron visualizer is included as `singleTetSubgrid` -- and another UI for a set of four tets sharing an edge in `ringTetSubgrid`.

## Quick examples

The executables `subgrid` and `dualSubgrid` accept three input types and write the
result to exactly the `-o` path you give. Run these from the repo root after
[building](#building):

```sh
# 1. Triangle mesh (OBJ / PLY / OFF), extracted on a 64³ grid
./build/subgrid -i ./data/meshes/spot.obj -r 64 -o ./out/spot.obj

# 2. Built-in signed distance function (see deps/sdf-dataset)
./build/subgrid -s Sphere -r 64 -o ./out/sphere.obj

# 3. Precomputed edge intersections (.npz): explicit tet mesh + hits, no grid
./build/subgrid --npz ./data/npz/wine_glass_N64_explicit.npz -o ./out/explicit.obj
```

A few small sample inputs ship in `data/` — `cube`, `spot`, and `wine_glass`,
each as both a mesh (`meshes/*.obj`) and a precomputed `.npz` (`npz/*_N64_explicit.npz`).

Swap `subgrid` → `dualSubgrid` for the dual (QEF) mesh from any of these inputs. 
See [Usage](#usage) for the full flag reference and [explicit input format](docs/explicit_input_format.md) for the
`.npz` layout.


## Dependencies

All dependencies are either fetched automatically by CMake (via CPM) 
or exist as submodules (`geometry-central`, `sdf-dataset`):

- [geometry-central](https://github.com/nmwsharp/geometry-central) — mesh data structures (submodule)
- [sdf-dataset](https://github.com/GeometryCollective/sdf-dataset) — built-in signed distance functions (submodule)
- [polyscope](https://github.com/nmwsharp/polyscope) — visualization (fetched by CMake; **optional**, see below)
- [fcpw](https://github.com/rohan-sawhney/fcpw) — BVH ray queries (fetched by CMake)
- [Eigen](https://eigen.tuxfamily.org) — linear algebra (via geometry-central)
- C++17 compiler

## Cloning

```sh
git clone --recursive git@github.com:hbaktash/subgrid-marching.git
cd subgrid-marching
git submodule update --init --recursive
```

## Building

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

All executables land in `build/`.

### Headless build (no Polyscope)

The `subgrid` and `dualSubgrid` CLIs write their output meshes without ever
opening a window. If you only need the CLI — e.g. on a server — build without the
viewer to drop the Polyscope / OpenGL / GLFW dependency entirely:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSUBGRID_POLYSCOPE_VIEWER=OFF
cmake --build build -j
```

This builds only `subgrid`, `dualSubgrid`, and the tests. The interactive demos
(`singleTetSubgrid`, `ringTetSubgrid`) and the on-screen viewer require the
default `-DSUBGRID_POLYSCOPE_VIEWER=ON`.

### Just want the algorithm? (no build, no dependencies)

If you only need the **per-tet** primal/dual local constructions and don't want
to build this project at all, use the copy-pasteable ports in
[`standalone/`](standalone/): a single header-only C++ file (STL only) or a
single pure-Python file. Drop one file into your project and call
`subgrid_primal` / `subgrid_dual`. See [`standalone/README.md`](standalone/README.md).

## Usage

Both `subgrid` and `dualSubgrid` take exactly one type of input: a precomputed `.npz` of edge
intersections (`--npz`), a triangle mesh file (`-i`), or a named built-in SDF (`-s`). 
With `--npz` the tet mesh and intersections come entirely from the file, so `-r` does not apply — see
[explicit input format](docs/explicit_input_format.md) for the format and
[`data/npz/`](data/npz/) for examples.

For what each pipeline accepts and produces per tet, see the
[construction policy](docs/construction_policy.md).

### subgrid — primal extraction

```
./build/subgrid -i <mesh.obj|ply|off>  -r <resolution>  -o <output.obj>  [options]
./build/subgrid -s <SDF name>          -r <resolution>  -o <output.obj>  [options]
./build/subgrid --npz <hits.npz>                          -o <output.obj>  [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | — | Input triangle mesh |
| `-s, --inputSDF` | — | Named built-in SDF, e.g. `Sphere`, `Torus` (ported from [sdf-dataset](https://github.com/GeometryCollective/sdf-dataset)) |
| `--npz` | — | Explicit tet mesh + precomputed edge intersections (`.npz`); `-r` is ignored |
| `-r, --tetGridResolution` | 64 | Grid resolution — builds an n³ cube grid (5n³ tetrahedra) |
| `-o, --output` | — | Output mesh path; if omitted, result is only visualized |

**Vertex merging** is set by two flags: the default is exact **combinatorial**
merge; `--mergeEPS <eps>` (with `eps ≥ 0`) switches to positional **numerical**
merge, and `--noMerge` emits the raw per-tet triangle soup. The trade-offs are in
[vertex merging](docs/combinatorial_merging.md).

An alternative **greedy** primal construction (fan/spiral triangulation of each
boundary curve) is available with `--greedy` -- intersection-free property is not
guaranteed in this mode.

Remaining options (`--mod2`, `--noViz`, `--noPBar`, `--inputSaveDir`) are listed
by `./build/subgrid -h`. The output mesh is written to exactly the `-o` path; the
Polyscope window opens after extraction unless `--noViz` is given.

### dualSubgrid — dual extraction

```
./build/dualSubgrid -i <mesh>  -r <resolution>  -o <output.obj>  [options]
./build/dualSubgrid -s <SDF>   -r <resolution>  -o <output.obj>  [options]
./build/dualSubgrid --npz <hits.npz>              -o <output.obj>  [options]
```

> **Note:** When using `--npz` input, the `.npz` file **must include per-intersection normal information** 
in addition to the standard edge intersection data. See [explicit input format documentation](docs/explicit_input_format.md) for details.


Flags shared with `subgrid` (`-i`, `-s`, `--npz`, `-r`, `-o`, `--mod2`, `--noViz`,
`--noPBar`, `--inputSaveDir`) behave identically.

Additional flags:

| Flag | Default | Description |
|------|---------|-------------|
| `-a, --alpha` | 0.1 | QEF regularization — higher values pull dual vertices toward the polygon centroid |
| `--pd, --projectDuals` | off | Clip each dual vertex back inside its local grid cell |

The dual mesh is written to exactly the `-o` path.

### singleTetSubgrid — interactive visualizer

```sh
./build/singleTetSubgrid [-o output_prefix.obj]
```

Opens a Polyscope GUI for building and visualizing normal surfaces inside a
single tetrahedron. Sliders control normal coordinates (corner and diagonal
cuts) and edge intersection counts directly.

## Input preprocessing

When a mesh file is provided, the input is automatically:
1. Welded (vertices with identical positions merged) and triangulated, so that
   "soup" representations of watertight meshes stay watertight.
2. Centered and scaled to fit within 98% of the unit cube `[-1, 1]³`.
3. Rigidly translated by a small fixed-seed offset (magnitude 0.01) to break
   axis-alignment with the grid and avoid degenerate ray/grid coincidences. The
   offset is half the grid-boundary clearance, so the mesh always stays inside
   the grid. Unlike per-vertex jitter, this keeps the mesh geometrically exact.

Use `--inputSaveDir <path>` to save and inspect the preprocessed mesh.

## A note on non-even (open) tets

`subgrid` and `dualSubgrid` report a `non-even tets` count. Each active tetrahedron should
see the isosurface cross its faces an even number of times (the "even-sum"
condition); when it holds, the per-tet boundary curves close up and the output is
watertight and orientable. A non-zero count means some tets had an odd crossing
parity, producing open curves — the output may then have small holes or be
non-orientable near those tets.

The dominant cause is a **non-watertight input**: meshes with holes/boundaries
will always report some non-even tets, and that is expected. A
watertight input (even if self-intersecting) — even one stored as a triangle soup — should report
`non-even tets: 0`. A small residual count on an otherwise-closed mesh usually
reflects near-degenerate ray/grid intersections at that resolution; trying a
different resolution typically changes it. A more robust ray intersection implementation 
in the future should resolve such cases.

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Baktash:2026:SMT,
  author    = {Baktash, Hossein and Gillespie, Mark and Crane, Keenan},
  title     = {Subgrid Marching Tetrahedra},
  journal   = {ACM Trans. Graph.},
  issue_date = {July 2026},
  volume    = {45},
  number    = {4},
  articleno = {57},
  numpages  = {20},
  year      = {2026},
  publisher = {Association for Computing Machinery},
  address   = {New York, NY, USA},
  doi       = {10.1145/3811358},
  url       = {https://doi.org/10.1145/3811358}
}
```

## License

This project is released under the [MIT License](LICENSE).

Dependencies are subject to their own licenses (all permissive — MIT, BSD, or
similar). See each dependency's repository for details.
