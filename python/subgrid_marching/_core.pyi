"""Type stubs for the compiled subgrid_marching._core extension."""

from __future__ import annotations

from typing import Any, Iterable, List, Sequence, Union

import numpy as np
import numpy.typing as npt

ArrayLike = Union[npt.NDArray[Any], Sequence[Any]]
EdgeTs = Union[npt.NDArray[np.float64], Sequence[Sequence[float]]]
MergeMode = str  # "combinatorial" | "numerical" | "none"

SIG_I: int
SIG_J: int
SIG_ORDER: int
SIG_EDGE_INT: int
SIG_K: int
SIG_TYPE: int

EDGE_VERTEX_PAIRS: npt.NDArray[np.int64]

class CombVertexSigType:
    NORMAL: CombVertexSigType
    SCOOP_FACE_STEINER: CombVertexSigType
    SCOOP_INTERIOR_STEINER: CombVertexSigType
    STANDALONE: CombVertexSigType
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...
    def __int__(self) -> int: ...

class Soup:
    @property
    def vertices(self) -> npt.NDArray[np.float64]: ...
    @property
    def face_offsets(self) -> npt.NDArray[np.int64]: ...
    @property
    def face_vertices(self) -> npt.NDArray[np.int64]: ...
    @property
    def signatures(self) -> npt.NDArray[np.int32]: ...
    @property
    def faces(self) -> List[npt.NDArray[np.int64]]: ...
    @property
    def num_faces(self) -> int: ...

class PrimalResult(Soup):
    @property
    def non_even(self) -> bool: ...

class DualResult(Soup):
    @property
    def dual_positions(self) -> npt.NDArray[np.float64]: ...
    @property
    def face_edge_tet_faces(self) -> npt.NDArray[np.int64]: ...
    @property
    def non_even(self) -> bool: ...
    @property
    def non_normal(self) -> bool: ...

class PipelineResult(Soup):
    @property
    def non_zero_tets(self) -> int: ...
    @property
    def non_normal_tets(self) -> int: ...
    @property
    def non_even_tets(self) -> int: ...
    @property
    def total_tets(self) -> int: ...
    @property
    def isect_time(self) -> float: ...
    @property
    def construction_time(self) -> float: ...
    @property
    def assembly_time(self) -> float: ...

class DualPipelineResult(PipelineResult):
    @property
    def assembled(self) -> bool: ...
    @property
    def dual_positions(self) -> npt.NDArray[np.float64]: ...
    @property
    def face_edge_tet_faces(self) -> npt.NDArray[np.int64]: ...

def edge_pair_to_index(i: int, j: int) -> int: ...
def available_sdfs() -> List[str]: ...

# ---- per-tet constructions -------------------------------------------------

def subgrid_primal(
    tet_positions: ArrayLike,
    tet_global_indices: ArrayLike,
    edge_isect_ts: EdgeTs,
    *,
    scoop_mid_vertices: bool = ...,
    scoop_bulge: float = ...,
) -> PrimalResult: ...
def subgrid_greedy(
    tet_positions: ArrayLike,
    tet_global_indices: ArrayLike,
    edge_isect_ts: EdgeTs,
    *,
    scoop_mid_vertices: bool = ...,
    scoop_bulge: float = ...,
) -> PrimalResult: ...
def subgrid_dual(
    tet_positions: ArrayLike,
    tet_global_indices: ArrayLike,
    edge_isect_ts: EdgeTs,
    edge_isect_normals: Union[Iterable[ArrayLike], None] = ...,
    *,
    reg_alpha: float = ...,
    project_duals: bool = ...,
    use_normals: bool = ...,
) -> DualResult: ...

# ---- primal pipelines ------------------------------------------------------

def primal_from_mesh(
    vertices: ArrayLike,
    faces: ArrayLike,
    resolution: int = ...,
    *,
    preprocess: bool = ...,
    seed: int = ...,
    cgal: bool = ...,
    mod2: bool = ...,
    greedy: bool = ...,
    scoop_bulge: float = ...,
    scoop_mid_vertices: bool = ...,
    merge: MergeMode = ...,
    merge_eps: float = ...,
    orient: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
    verbose: bool = ...,
) -> PipelineResult: ...
def primal_from_mesh_file(
    path: str,
    resolution: int = ...,
    *,
    preprocess: bool = ...,
    seed: int = ...,
    cgal: bool = ...,
    mod2: bool = ...,
    greedy: bool = ...,
    scoop_bulge: float = ...,
    scoop_mid_vertices: bool = ...,
    merge: MergeMode = ...,
    merge_eps: float = ...,
    orient: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
    verbose: bool = ...,
) -> PipelineResult: ...
def primal_from_sdf(
    name: str,
    resolution: int = ...,
    *,
    step_size: float = ...,
    mod2: bool = ...,
    greedy: bool = ...,
    scoop_bulge: float = ...,
    scoop_mid_vertices: bool = ...,
    merge: MergeMode = ...,
    merge_eps: float = ...,
    orient: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
    verbose: bool = ...,
) -> PipelineResult: ...
def primal_from_npz(
    path: str,
    *,
    mod2: bool = ...,
    greedy: bool = ...,
    scoop_bulge: float = ...,
    scoop_mid_vertices: bool = ...,
    merge: MergeMode = ...,
    merge_eps: float = ...,
    orient: bool = ...,
    progress: bool = ...,
    verbose: bool = ...,
) -> PipelineResult: ...

# ---- dual pipelines --------------------------------------------------------

def dual_from_mesh(
    vertices: ArrayLike,
    faces: ArrayLike,
    resolution: int = ...,
    *,
    preprocess: bool = ...,
    seed: int = ...,
    cgal: bool = ...,
    mod2: bool = ...,
    reg_alpha: float = ...,
    project_duals: bool = ...,
    use_normals: bool = ...,
    assemble: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
) -> DualPipelineResult: ...
def dual_from_mesh_file(
    path: str,
    resolution: int = ...,
    *,
    preprocess: bool = ...,
    seed: int = ...,
    cgal: bool = ...,
    mod2: bool = ...,
    reg_alpha: float = ...,
    project_duals: bool = ...,
    use_normals: bool = ...,
    assemble: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
) -> DualPipelineResult: ...
def dual_from_sdf(
    name: str,
    resolution: int = ...,
    *,
    step_size: float = ...,
    mod2: bool = ...,
    reg_alpha: float = ...,
    project_duals: bool = ...,
    use_normals: bool = ...,
    assemble: bool = ...,
    robust: bool = ...,
    progress: bool = ...,
) -> DualPipelineResult: ...
def dual_from_npz(
    path: str,
    *,
    mod2: bool = ...,
    reg_alpha: float = ...,
    project_duals: bool = ...,
    use_normals: bool = ...,
    assemble: bool = ...,
    progress: bool = ...,
) -> DualPipelineResult: ...
