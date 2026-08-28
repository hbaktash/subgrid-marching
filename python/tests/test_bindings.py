"""Tests for the subgrid_marching Python bindings.

Run from a repo checkout after `pip install .`:

    pytest python/tests

The parity tests compare the bindings against standalone/python/subgrid_mt.py — the
dependency-free reference port of the same per-tet constructions — and the pipeline
tests use the sample inputs in data/. Each group skips itself if what it needs is
not present, so the suite also runs against a bare wheel.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import subgrid_marching as smt

REPO_ROOT = Path(__file__).resolve().parents[2]
PORT_PATH = REPO_ROOT / "standalone" / "python" / "subgrid_mt.py"
MESH_DIR = REPO_ROOT / "data" / "meshes"
NPZ_DIR = REPO_ROOT / "data" / "npz"

UNIT_TET = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
IDENTITY_INDICES = np.array([0, 1, 2, 3])

# One t per edge (0,1), (0,2), (0,3): the surface cuts off the corner at vertex 0.
CORNER_V0 = [[0.3], [0.4], [0.5], [], [], []]
# Two hits on each of the three edges at vertex 0: a "double corner".
DOUBLE_V0 = [[0.2, 0.6], [0.3, 0.7], [0.4, 0.8], [], [], []]
# A quad separating {0,1} from {2,3}.
QUAD_01_23 = [[], [0.5], [0.5], [0.5], [0.5], []]


@pytest.fixture(scope="module")
def port():
    """The dependency-free reference port, imported straight from the repo."""
    if not PORT_PATH.exists():
        pytest.skip(f"reference port not found at {PORT_PATH}")
    spec = importlib.util.spec_from_file_location("subgrid_mt_reference", PORT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def as_xyz(points):
    """The port returns Vec3 objects; the bindings return an (N, 3) array."""
    return np.array([[p.x, p.y, p.z] for p in points], dtype=float).reshape(-1, 3)


def face_list(result):
    return [list(f) for f in result.faces]


# ---------------------------------------------------------------------------
# conventions
# ---------------------------------------------------------------------------


def test_edge_order_is_the_documented_one():
    expected = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    assert [tuple(p) for p in smt.EDGE_VERTEX_PAIRS] == expected
    for index, (i, j) in enumerate(expected):
        assert smt.edge_pair_to_index(i, j) == index
        assert smt.edge_pair_to_index(j, i) == index


def test_signature_columns_are_distinct():
    columns = [smt.SIG_I, smt.SIG_J, smt.SIG_ORDER, smt.SIG_EDGE_INT, smt.SIG_K, smt.SIG_TYPE]
    assert columns == [0, 1, 2, 3, 4, 5]


def test_available_sdfs_is_non_empty():
    names = smt.available_sdfs()
    assert names == sorted(names)
    assert "Sphere" in names


# ---------------------------------------------------------------------------
# per-tet primal
# ---------------------------------------------------------------------------


def test_corner_cut_is_one_triangle_through_the_hit_points():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, CORNER_V0)

    assert not result.non_even
    assert result.num_faces == 1
    assert result.vertices.shape == (3, 3)
    np.testing.assert_allclose(
        np.sort(result.vertices, axis=0),
        np.sort(np.array([[0.3, 0, 0], [0, 0.4, 0], [0, 0, 0.5]]), axis=0),
    )
    assert sorted(face_list(result)[0]) == [0, 1, 2]


def test_empty_tet_gives_an_empty_soup():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, [[], [], [], [], [], []])

    assert result.vertices.shape == (0, 3)
    assert result.num_faces == 0
    assert result.faces == []
    assert result.face_offsets.shape in {(0,), (1,)}


def test_quad_configuration_produces_a_closed_ring():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, QUAD_01_23)

    assert not result.non_even
    assert result.vertices.shape[0] >= 4
    assert result.num_faces >= 1


def test_odd_edge_parity_is_reported_as_non_even():
    # A single intersection on one edge cannot close up into a loop.
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, [[0.5], [], [], [], [], []])
    assert result.non_even


def test_greedy_matches_the_even_sum_construction_on_a_corner():
    even_sum = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, CORNER_V0)
    greedy = smt.subgrid_greedy(UNIT_TET, IDENTITY_INDICES, CORNER_V0)

    np.testing.assert_allclose(np.sort(greedy.vertices, axis=0), np.sort(even_sum.vertices, axis=0))
    assert greedy.num_faces == even_sum.num_faces


# ---------------------------------------------------------------------------
# array plumbing
# ---------------------------------------------------------------------------


def test_arrays_have_the_documented_dtypes_and_shapes():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, DOUBLE_V0)

    assert result.vertices.dtype == np.float64
    assert result.vertices.ndim == 2 and result.vertices.shape[1] == 3
    assert result.face_offsets.dtype == np.int64
    assert result.face_vertices.dtype == np.int64
    assert result.signatures.dtype == np.int32
    assert result.signatures.shape == (result.vertices.shape[0], 6)


def test_returned_arrays_own_their_buffer_rather_than_copying():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, DOUBLE_V0)

    # `base` is the capsule holding the moved C++ vector: the data was handed over,
    # not copied out.
    assert result.vertices.base is not None
    assert result.signatures.base is not None


def test_faces_views_agree_with_the_csr_arrays():
    result = smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, DOUBLE_V0)

    offsets = result.face_offsets
    assert offsets[0] == 0
    assert offsets[-1] == result.face_vertices.shape[0]
    assert len(result.faces) == result.num_faces == offsets.shape[0] - 1

    for f, face in enumerate(result.faces):
        np.testing.assert_array_equal(face, result.face_vertices[offsets[f] : offsets[f + 1]])
        assert len(face) >= 3
        # The views share the CSR buffer instead of copying it.
        assert face.base is result.face_vertices

    assert result.face_vertices.max() < result.vertices.shape[0]


def test_signatures_identify_the_grid_edge_a_vertex_sits_on():
    globals_ = np.array([7, 11, 13, 17])
    result = smt.subgrid_primal(UNIT_TET, globals_, CORNER_V0)

    # Every vertex of a corner cut lies on an edge out of local vertex 0 -> global 7,
    # is the first hit on that edge, and is an ordinary (mergeable) vertex.
    for signature in result.signatures:
        assert signature[smt.SIG_I] == 7
        assert signature[smt.SIG_J] in (11, 13, 17)
        assert signature[smt.SIG_ORDER] == 1
        assert signature[smt.SIG_TYPE] == int(smt.CombVertexSigType.NORMAL)


# ---------------------------------------------------------------------------
# per-tet dual
# ---------------------------------------------------------------------------


def test_dual_returns_one_dual_point_per_boundary_loop():
    normals = [np.tile([0.0, 0.0, 1.0], (len(ts), 1)) for ts in CORNER_V0]
    result = smt.subgrid_dual(UNIT_TET, IDENTITY_INDICES, CORNER_V0, normals)

    assert result.num_faces == 1
    assert result.dual_positions.shape == (1, 3)
    assert result.face_edge_tet_faces.shape == (result.face_vertices.shape[0], 3)
    assert not result.non_even and not result.non_normal


def test_dual_without_normals_uses_the_boundary_centroid():
    result = smt.subgrid_dual(UNIT_TET, IDENTITY_INDICES, CORNER_V0, use_normals=False)

    np.testing.assert_allclose(result.dual_positions[0], result.vertices.mean(axis=0))


def test_dual_face_tags_name_tet_faces_of_the_loop_edges():
    result = smt.subgrid_dual(UNIT_TET, IDENTITY_INDICES, CORNER_V0, use_normals=False)

    tags = result.face_edge_tet_faces
    assert tags.shape[0] == result.face_vertices.shape[0]
    for tag in tags:
        assert len(set(tag.tolist())) == 3
        assert set(tag.tolist()) <= {0, 1, 2, 3}


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------


def test_rejects_a_badly_shaped_tet():
    with pytest.raises(ValueError, match=r"\(4, 3\)"):
        smt.subgrid_primal(np.zeros((3, 3)), IDENTITY_INDICES, CORNER_V0)


def test_rejects_the_wrong_number_of_edges():
    with pytest.raises(ValueError, match="6 entries"):
        smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, [[0.5]] * 5)


def test_rejects_unsorted_t_values():
    with pytest.raises(ValueError, match="sorted ascending"):
        smt.subgrid_primal(UNIT_TET, IDENTITY_INDICES, [[0.6, 0.2], [], [], [], [], []])


def test_rejects_normals_that_do_not_match_the_intersections():
    with pytest.raises(ValueError, match="normals"):
        smt.subgrid_dual(UNIT_TET, IDENTITY_INDICES, CORNER_V0, [np.zeros((2, 3))] * 6)


def test_dual_requires_normals_unless_they_are_switched_off():
    with pytest.raises(ValueError, match="use_normals=False"):
        smt.subgrid_dual(UNIT_TET, IDENTITY_INDICES, CORNER_V0)


def test_rejects_an_unknown_merge_mode():
    pytest.importorskip("numpy")
    with pytest.raises(ValueError, match="combinatorial"):
        smt.primal_from_mesh(np.zeros((3, 3)), np.array([[0, 1, 2]]), 4, merge="sometimes")


# ---------------------------------------------------------------------------
# parity with the dependency-free reference port
# ---------------------------------------------------------------------------

PARITY_BULGE = 1e-4  # standalone/'s default; the bindings default to the pipeline's 1e-3


def random_configs(seed, count, max_hits=3):
    rng = np.random.default_rng(seed)
    for _ in range(count):
        counts = rng.integers(0, max_hits + 1, size=6)
        ts = [sorted(rng.random(int(c)).tolist()) for c in counts]
        base = int(rng.integers(0, 50))
        indices = rng.permutation(np.arange(4) + base).tolist()
        yield ts, indices, rng


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_primal_matches_the_reference_port(port, seed):
    for ts, indices, _ in random_configs(seed, 40):
        mine = smt.subgrid_primal(UNIT_TET, np.array(indices), ts, scoop_bulge=PARITY_BULGE)
        theirs = port.subgrid_primal(
            [tuple(p) for p in UNIT_TET], indices, ts, scoop_mid_vertex_bulge=PARITY_BULGE
        )

        np.testing.assert_allclose(mine.vertices, as_xyz(theirs.vertices), atol=1e-12)
        assert face_list(mine) == [list(f) for f in theirs.faces]
        assert mine.non_even == theirs.non_even

        theirs_signatures = np.array(
            [[s.i, s.j, s.order, s.k, int(s.type)] for s in theirs.signatures], dtype=np.int32
        ).reshape(-1, 5)
        np.testing.assert_array_equal(
            mine.signatures[:, [smt.SIG_I, smt.SIG_J, smt.SIG_ORDER, smt.SIG_K, smt.SIG_TYPE]],
            theirs_signatures,
        )


@pytest.mark.parametrize("seed", [3, 4])
def test_dual_matches_the_reference_port(port, seed):
    for ts, indices, rng in random_configs(seed, 30):
        normals = []
        for edge_ts in ts:
            n = rng.normal(size=(len(edge_ts), 3))
            if len(edge_ts):
                n /= np.linalg.norm(n, axis=1, keepdims=True)
            normals.append(n)

        mine = smt.subgrid_dual(UNIT_TET, np.array(indices), ts, normals)
        theirs = port.subgrid_dual(
            [tuple(p) for p in UNIT_TET],
            indices,
            ts,
            [[tuple(v) for v in edge] for edge in normals],
        )

        np.testing.assert_allclose(mine.vertices, as_xyz(theirs.vertices), atol=1e-12)
        assert face_list(mine) == [list(f) for f in theirs.faces]
        # The port hand-rolls the 3x3 QEF solve instead of using Eigen.
        np.testing.assert_allclose(mine.dual_positions, as_xyz(theirs.dual_positions), atol=1e-6)
        assert mine.non_even == theirs.non_even
        assert mine.non_normal == theirs.non_normal


# ---------------------------------------------------------------------------
# pipelines
# ---------------------------------------------------------------------------


def polygon_edges(result):
    """Undirected (min, max) vertex pairs of every polygon edge in a CSR soup."""
    offsets = result.face_offsets
    edges = []
    for f in range(result.num_faces):
        loop = result.face_vertices[offsets[f] : offsets[f + 1]]
        edges.extend(
            (min(int(a), int(b)), max(int(a), int(b)))
            for a, b in zip(loop, np.roll(loop, -1))
        )
    return edges


def read_obj(path):
    vertices, faces = [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "v":
            vertices.append([float(x) for x in parts[1:4]])
        elif parts[0] == "f":
            faces.append([int(x.split("/")[0]) - 1 for x in parts[1:]])
    return np.array(vertices), np.array(faces)


@pytest.fixture(scope="module")
def spot():
    path = MESH_DIR / "spot.obj"
    if not path.exists():
        pytest.skip(f"sample mesh not found at {path}")
    return path


def test_primal_pipeline_on_a_watertight_mesh(spot):
    result = smt.primal_from_mesh_file(str(spot), resolution=16)

    assert result.non_even_tets == 0, "a watertight input must satisfy the even-sum condition"
    assert 0 < result.non_zero_tets < result.total_tets
    assert result.vertices.shape[0] > 0 and result.num_faces > 0
    assert result.signatures.shape == (result.vertices.shape[0], 6)

    # The output of a watertight input must itself be closed: every edge is shared
    # by exactly two faces. This is what combinatorial merging buys.
    counts = Counter(polygon_edges(result))
    assert set(counts.values()) == {2}


def test_in_memory_mesh_matches_the_file_path(spot):
    vertices, faces = read_obj(spot)

    from_arrays = smt.primal_from_mesh(vertices, faces, resolution=16)
    from_file = smt.primal_from_mesh_file(str(spot), resolution=16)

    np.testing.assert_array_equal(from_arrays.vertices, from_file.vertices)
    np.testing.assert_array_equal(from_arrays.face_vertices, from_file.face_vertices)


def test_merge_modes_change_only_the_vertex_count(spot):
    merged = smt.primal_from_mesh_file(str(spot), resolution=16)
    raw = smt.primal_from_mesh_file(str(spot), resolution=16, merge="none")

    assert raw.vertices.shape[0] > merged.vertices.shape[0]
    assert raw.num_faces == merged.num_faces
    assert raw.signatures.shape[0] == 0  # no merging, so no signatures


def test_dual_pipeline_assembled_and_raw(spot):
    assembled = smt.dual_from_mesh_file(str(spot), resolution=16)
    assert assembled.assembled
    assert assembled.non_even_tets == 0
    assert assembled.vertices.shape[0] > 0 and assembled.num_faces > 0

    raw = smt.dual_from_mesh_file(str(spot), resolution=16, assemble=False)
    assert not raw.assembled
    # One dual point per boundary loop, one tet-face tag per polygon edge.
    assert raw.dual_positions.shape == (raw.num_faces, 3)
    assert raw.face_edge_tet_faces.shape == (raw.face_vertices.shape[0], 3)
    # The assembled dual has one vertex per loop of the raw soup.
    assert assembled.vertices.shape[0] == raw.num_faces


def test_sdf_pipeline():
    result = smt.primal_from_sdf("Sphere", resolution=16)

    assert result.non_even_tets == 0
    assert result.vertices.shape[0] > 0
    radii = np.linalg.norm(result.vertices, axis=1)
    assert radii.min() > 0.1  # a sphere, not a degenerate blob


def test_npz_pipeline():
    path = NPZ_DIR / "wine_glass_N64_explicit.npz"
    if not path.exists():
        pytest.skip(f"sample npz not found at {path}")

    result = smt.primal_from_npz(str(path))
    assert result.vertices.shape[0] > 0 and result.num_faces > 0

    # This archive carries no isect_normals, so the dual needs use_normals=False.
    with pytest.raises(RuntimeError, match="use_normals=False"):
        smt.dual_from_npz(str(path))
    centroid_dual = smt.dual_from_npz(str(path), use_normals=False)
    assert centroid_dual.vertices.shape[0] > 0
