"""Subgrid Marching Tetrahedra — intersection-free isosurface extraction.

Implementation of `Subgrid Marching Tetrahedra
<https://hbaktash.github.io/projects/subgrid-marching-tetrahedra/index.html>`_,
which extracts isosurfaces from the intersections between a surface and the edges
of a tetrahedral grid. The output is guaranteed free of self-intersections.

Two layers are exposed, both backed by the C++ core in ``src/subgrid_MT/``:

Per tetrahedron
    :func:`subgrid_primal`, :func:`subgrid_dual` and :func:`subgrid_greedy` build the
    local surface inside one tet from its edge intersections. Use these when you have
    your own tet mesh and intersection data; merge across tets with the per-vertex
    ``signatures`` (see ``examples/global_merge.py``).

Whole pipelines
    ``primal_from_mesh`` / ``primal_from_mesh_file`` / ``primal_from_sdf`` /
    ``primal_from_npz`` and their ``dual_from_*`` counterparts run the full extraction
    over a grid, mirroring the ``subgrid`` and ``dualSubgrid`` command line tools.

Quick start::

    import subgrid_marching as smt

    result = smt.primal_from_mesh_file("data/meshes/spot.obj", resolution=64)
    result.vertices       # (N, 3) float64
    result.faces          # list of index arrays, one per face
    result.non_even_tets  # 0 on a watertight input
"""

from __future__ import annotations

from ._core import (
    EDGE_VERTEX_PAIRS,
    SIG_EDGE_INT,
    SIG_I,
    SIG_J,
    SIG_K,
    SIG_ORDER,
    SIG_TYPE,
    CombVertexSigType,
    DualPipelineResult,
    DualResult,
    PipelineResult,
    PrimalResult,
    Soup,
    available_sdfs,
    dual_from_mesh,
    dual_from_mesh_file,
    dual_from_npz,
    dual_from_sdf,
    edge_pair_to_index,
    primal_from_mesh,
    primal_from_mesh_file,
    primal_from_npz,
    primal_from_sdf,
    subgrid_dual,
    subgrid_greedy,
    subgrid_primal,
)

try:  # populated by pip; absent when importing straight out of a build tree
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version("subgrid-marching")
    except PackageNotFoundError:  # pragma: no cover
        __version__ = "0.0.0+local"
except ImportError:  # pragma: no cover - Python < 3.8
    __version__ = "0.0.0+local"

__all__ = [
    # per-tet constructions
    "subgrid_primal",
    "subgrid_dual",
    "subgrid_greedy",
    # pipelines
    "primal_from_mesh",
    "primal_from_mesh_file",
    "primal_from_sdf",
    "primal_from_npz",
    "dual_from_mesh",
    "dual_from_mesh_file",
    "dual_from_sdf",
    "dual_from_npz",
    "available_sdfs",
    # result types
    "Soup",
    "PrimalResult",
    "DualResult",
    "PipelineResult",
    "DualPipelineResult",
    # conventions
    "CombVertexSigType",
    "EDGE_VERTEX_PAIRS",
    "edge_pair_to_index",
    "SIG_I",
    "SIG_J",
    "SIG_ORDER",
    "SIG_EDGE_INT",
    "SIG_K",
    "SIG_TYPE",
    "__version__",
]
