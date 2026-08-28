"""Extract an isosurface from a mesh and write it out as an OBJ.

Mirrors `subgrid -i <mesh> -r <resolution> -o <out.obj>`. Run from the repo root:

    python python/examples/extract_and_save.py data/meshes/spot.obj out/spot_python.obj
"""

import sys
from pathlib import Path

import numpy as np

import subgrid_marching as smt


def write_obj(path, vertices, face_offsets, face_vertices):
    """Minimal OBJ writer; faces may be n-gons, so walk the CSR arrays."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for x, y, z in vertices:
            out.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for f in range(len(face_offsets) - 1):
            loop = face_vertices[face_offsets[f] : face_offsets[f + 1]]
            out.write("f " + " ".join(str(int(i) + 1) for i in loop) + "\n")


def main():
    mesh_path = sys.argv[1] if len(sys.argv) > 1 else "data/meshes/spot.obj"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "out/subgrid_python.obj"
    resolution = int(sys.argv[3]) if len(sys.argv) > 3 else 64

    result = smt.primal_from_mesh_file(mesh_path, resolution=resolution, progress=True)

    print(f" non-zero tets:   {result.non_zero_tets}/{result.total_tets}")
    print(f" non-even tets:   {result.non_even_tets}/{result.non_zero_tets}")
    if result.non_even_tets:
        print(
            " warning: the even-sum condition failed somewhere, so the output may have "
            "holes or pinch vertices. Usually a non-watertight input; try another "
            "resolution or seed."
        )
    print(f" output: {result.vertices.shape[0]} verts, {result.num_faces} faces")
    print(f" query: {result.isect_time:.3f}s  construction: {result.construction_time:.3f}s")

    write_obj(out_path, result.vertices, result.face_offsets, result.face_vertices)
    print(f" saved to {out_path}")

    # The dual (QEF) surface from the same input is one call away.
    dual = smt.dual_from_mesh_file(mesh_path, resolution=resolution)
    print(f" dual: {dual.vertices.shape[0]} verts, {dual.num_faces} faces")


if __name__ == "__main__":
    main()
