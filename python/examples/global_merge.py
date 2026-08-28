"""Stitching per-tet soups into one mesh, using the combinatorial signatures.

Each output vertex carries a signature that says which global grid edge it lies on and
where along that edge it sits. Two vertices produced by different tets are the same
point exactly when their signatures match — no epsilon, no positional search. The one
exception is the STANDALONE type, which must never merge: give each a fresh index.

    python python/examples/global_merge.py
"""

import numpy as np

import subgrid_marching as smt

IDENTITY_COLUMNS = [smt.SIG_I, smt.SIG_J, smt.SIG_ORDER, smt.SIG_K, smt.SIG_TYPE]
STANDALONE = int(smt.CombVertexSigType.STANDALONE)


def merge_soups(results):
    """Merge per-tet results into one indexed mesh. Returns (vertices, faces)."""
    vertices = []
    faces = []
    by_signature = {}

    for result in results:
        signatures = result.signatures
        remap = np.empty(result.vertices.shape[0], dtype=np.int64)

        for local, position in enumerate(result.vertices):
            identity = tuple(signatures[local, IDENTITY_COLUMNS].tolist())
            standalone = signatures[local, smt.SIG_TYPE] == STANDALONE
            existing = None if standalone else by_signature.get(identity)

            if existing is None:
                remap[local] = len(vertices)
                if not standalone:
                    by_signature[identity] = remap[local]
                vertices.append(position)
            else:
                remap[local] = existing

        faces.extend(remap[face].tolist() for face in result.faces)

    return np.array(vertices).reshape(-1, 3), faces


# Two tets sharing the face (1, 2, 3), with the surface cutting the corner at global
# vertex 1. The shared edges (1,2) and (1,3) carry the same t in both tets, as they
# must: the intersections are a property of the edge, not of the tet you look from.
positions = np.array(
    [
        [0.0, 0.0, 0.0],  # 0
        [1.0, 0.0, 0.0],  # 1
        [0.0, 1.0, 0.0],  # 2
        [0.0, 0.0, 1.0],  # 3
        [1.0, 1.0, 1.0],  # 4
    ]
)

tets = [
    # (global vertex ids, per-edge t values in the canonical edge order)
    ([0, 1, 2, 3], [[0.7], [], [], [0.3], [0.3], []]),
    ([1, 2, 3, 4], [[0.3], [0.3], [0.3], [], [], []]),
]

results = [
    smt.subgrid_primal(positions[list(indices)], np.array(indices), edge_ts)
    for indices, edge_ts in tets
]

per_tet_vertices = sum(r.vertices.shape[0] for r in results)
vertices, faces = merge_soups(results)

print(f"per-tet vertices: {per_tet_vertices}")
print(f"merged vertices:  {vertices.shape[0]}  (the two points on the shared face are one)")
print(f"faces:            {faces}")
print(vertices)

# The merge is exact rather than approximate: the duplicates it removed were bitwise
# identical, because both tets interpolate the same t along the same global edge.
assert vertices.shape[0] == per_tet_vertices - 2
