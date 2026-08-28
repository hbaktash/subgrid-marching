"""One tetrahedron, one surface patch.

The per-tet entry points are the heart of the method: given where the surface crosses
the six edges of a tet, they build the local piece of the isosurface. Everything else
(the pipelines, the CLIs) is a loop over tets around these calls.

    python python/examples/single_tet.py
"""

import numpy as np

import subgrid_marching as smt

# The tet: four positions, and the ids these vertices have in your global tet mesh.
# The ids matter — they are what makes neighbouring tets agree on shared points.
tet_positions = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)
tet_global_indices = np.array([0, 1, 2, 3])

# Where the surface crosses each edge, in the canonical edge order
#   (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)     — smt.EDGE_VERTEX_PAIRS
# For edge (i, j) with i < j, t denotes (1 - t) * v_i + t * v_j, sorted ascending.
# Here one crossing on each edge out of vertex 0: the surface cuts off that corner.
edge_isect_ts = [[0.3], [0.4], [0.5], [], [], []]

primal = smt.subgrid_primal(tet_positions, tet_global_indices, edge_isect_ts)

print("primal:", primal)
print("  vertices:\n", primal.vertices)
print("  faces:", [face.tolist() for face in primal.faces])
print("  even-sum:", not primal.non_even)

# Faces are general n-gons, so they arrive in CSR form. `.faces` is a list of NumPy
# views over the same buffer; the raw arrays are there when you want them vectorized.
print("  face sizes:", np.diff(primal.face_offsets))

# Each vertex carries the combinatorial signature that identifies it across tets:
# [i, j, order, edge_int, k, type] — see examples/global_merge.py.
print("  signatures:\n", primal.signatures)

# The dual construction keeps the boundary loop instead of filling it, and places one
# point per loop by a QEF solve over the intersection points and their surface normals.
edge_isect_normals = [
    np.array([[1.0, 0.0, 0.0]]),
    np.array([[0.0, 1.0, 0.0]]),
    np.array([[0.0, 0.0, 1.0]]),
    np.zeros((0, 3)),
    np.zeros((0, 3)),
    np.zeros((0, 3)),
]

dual = smt.subgrid_dual(tet_positions, tet_global_indices, edge_isect_ts, edge_isect_normals)

print("\ndual:", dual)
print("  boundary loop:", [face.tolist() for face in dual.faces])
print("  dual point:", dual.dual_positions)

# Without normals the dual point falls back to the centroid of its boundary polygon.
centroid = smt.subgrid_dual(tet_positions, tet_global_indices, edge_isect_ts, use_normals=False)
print("  centroid instead of QEF:", centroid.dual_positions)
