"""subgrid_mt.py - self-contained, dependency-free port of the Subgrid Marching
Tetrahedra per-tet reconstruction (primal + dual).

Pure Python, standard library only (numpy is NOT required). Single file,
copy-pastable. This is a FAITHFUL transcription of the repository core in
src/subgrid_MT/ + src/nc/ + the pieces of common/utils it depends on. Do NOT
change the algorithm here; if the core changes, re-port from the maintained
core in src/subgrid_MT/.

Two entry points, matching docs/construction_policy.md "Input":
    subgrid_primal(tet_positions, tet_global_indices, edge_isect_ts)
    subgrid_dual  (tet_positions, tet_global_indices, edge_isect_ts,
                   edge_isect_normals, reg_alpha, project_duals)

Both return a per-tet soup with global-remapped combinatorial signatures (already
within-tet merged, as in the core). Global cross-tet merging is left to the
caller; use signature_key() to dedup (see standalone/examples).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import IntEnum
from typing import List, Optional, Tuple

# ============================================================================
# Vector3  (minimal port of geometrycentral::Vector3)
# ============================================================================


class Vec3:
    __slots__ = ("x", "y", "z")

    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @staticmethod
    def zero() -> "Vec3":
        return Vec3(0.0, 0.0, 0.0)

    @staticmethod
    def undefined() -> "Vec3":
        n = float("nan")
        return Vec3(n, n, n)

    def is_defined(self) -> bool:
        return not (math.isnan(self.x) or math.isnan(self.y) or math.isnan(self.z))

    def __add__(self, o: "Vec3") -> "Vec3":
        return Vec3(self.x + o.x, self.y + o.y, self.z + o.z)

    def __sub__(self, o: "Vec3") -> "Vec3":
        return Vec3(self.x - o.x, self.y - o.y, self.z - o.z)

    def __mul__(self, s: float) -> "Vec3":
        return Vec3(self.x * s, self.y * s, self.z * s)

    __rmul__ = __mul__

    def __truediv__(self, s: float) -> "Vec3":
        return Vec3(self.x / s, self.y / s, self.z / s)

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z

    def __getitem__(self, i: int) -> float:
        return (self.x, self.y, self.z)[i]

    def norm2(self) -> float:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def norm(self) -> float:
        return math.sqrt(self.norm2())

    def __repr__(self) -> str:
        return f"Vec3({self.x}, {self.y}, {self.z})"


def dot(a: Vec3, b: Vec3) -> float:
    return a.x * b.x + a.y * b.y + a.z * b.z


def cross(a: Vec3, b: Vec3) -> Vec3:
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)


def _as_vec3(p) -> Vec3:
    if isinstance(p, Vec3):
        return p
    return Vec3(p[0], p[1], p[2])


# C++-style truncation-toward-zero integer division by 2 (matches int math in core).
def _tdiv2(a: int) -> int:
    return int(a / 2)


# ============================================================================
# tet constants + generic set/geometry helpers  (from common/utils)
# ============================================================================
ALL_TET_PAIRS: List[Tuple[int, int]] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
ALL_TET_TRIPLETS: List[Tuple[int, int, int]] = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]


def edge_pair_to_index(i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    table = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}
    if (i, j) in table:
        return table[(i, j)]
    raise ValueError("edge_pair_to_index: Invalid edge pair")


def complement(inputs) -> List[int]:
    universe = {0, 1, 2, 3}
    for x in inputs:
        universe.discard(x)
    return sorted(universe)


def all_pairs(universe) -> List[Tuple[int, int]]:
    u = list(universe)
    return [(u[a], u[b]) for a in range(len(u)) for b in range(a + 1, len(u))]


def triangle_inequality(a: int, b: int, c: int) -> Tuple[bool, int, int, int]:
    mx = max(a, b, c)
    feasible = (a + b + c) >= 2 * mx
    res_a = res_b = res_c = 0
    if not feasible:
        if a == mx:
            res_a = a - b - c
        elif b == mx:
            res_b = b - a - c
        else:
            res_c = c - a - b
    return feasible, res_a, res_b, res_c


def segment_triangle_intersection(p0: Vec3, p1: Vec3, a: Vec3, b: Vec3, c: Vec3):
    EPS = 1e-12
    d = p1 - p0
    edge1 = b - a
    edge2 = c - a
    pvec = cross(d, edge2)
    det = dot(edge1, pvec)
    if abs(det) < EPS:
        return False, 0.0
    inv_det = 1.0 / det
    tvec = p0 - a
    u = dot(tvec, pvec) * inv_det
    if u < -EPS or u > 1.0 + EPS:
        return False, 0.0
    qvec = cross(tvec, edge1)
    v = dot(d, qvec) * inv_det
    if v < -EPS or u + v > 1.0 + EPS:
        return False, 0.0
    t_ray = dot(edge2, qvec) * inv_det
    if d.norm2() <= EPS:
        return False, 0.0
    t_param = t_ray
    if t_param < -EPS or t_param > 1.0 + EPS:
        return False, 0.0
    return True, max(0.0, min(1.0, t_param))


def segment_tet_intersection(p0: Vec3, p1: Vec3, tet_positions: List[Vec3]):
    out_t = math.inf
    for triplet in ALL_TET_TRIPLETS:
        a, b, c = tet_positions[triplet[0]], tet_positions[triplet[1]], tet_positions[triplet[2]]
        hit, t_face = segment_triangle_intersection(p0, p1, a, b, c)
        if hit:
            out_t = min(out_t, t_face)
    return (out_t != math.inf), out_t


def is_inside_tet(p: Vec3, tet_positions: List[Vec3], eps: float = 1e-12) -> bool:
    def signed_tet_vol(a, b, c, d):
        return dot(cross(b - a, c - a), d - a)

    from_pd_vol = 0.0
    for triplet in ALL_TET_TRIPLETS:
        from_pd_vol += abs(signed_tet_vol(tet_positions[triplet[0]], tet_positions[triplet[1]],
                                          tet_positions[triplet[2]], p))
    tet_vol = abs(signed_tet_vol(tet_positions[0], tet_positions[1], tet_positions[2], tet_positions[3]))
    return from_pd_vol <= tet_vol + eps


def generate_uniform_isect_ts_single_edge(edge_int: int, t_max: float = 1.0) -> List[float]:
    return [float(i + 1) * t_max / float(edge_int + 1) for i in range(edge_int)]


# ============================================================================
# Normal-coordinate containers  (from nc/normalCoord_containers.h)
# ============================================================================
def _epair(i: int, j: int) -> Tuple[int, int]:
    return (i, j) if i <= j else (j, i)


class EdgeInts:
    def __init__(self, source=None):
        self.data = {}
        if source is None:
            return
        if isinstance(source, EdgeInts):
            self.data = dict(source.data)
        elif len(source) == 6 and all(isinstance(x, (int, float)) for x in source):
            # array<int,6> counts in canonical edge order
            for e, (i, j) in enumerate(ALL_TET_PAIRS):
                self.set(i, j, int(source[e]))
        elif len(source) == 6:
            # array of 6 lists (edge_isect_ts): use lengths
            for e, (i, j) in enumerate(ALL_TET_PAIRS):
                self.set(i, j, len(source[e]))
        else:
            raise ValueError("EdgeInts: unsupported source")

    def get(self, i: int, j: int) -> int:
        return self.data.get(_epair(i, j), 0)

    def set(self, i: int, j: int, v: int):
        self.data[_epair(i, j)] = int(v)

    def __getitem__(self, ij) -> int:
        return self.get(ij[0], ij[1])

    def inc(self, i: int, j: int):
        k = _epair(i, j)
        self.data[k] = self.data.get(k, 0) + 1

    def to_array(self) -> List[int]:
        return [self.get(i, j) for (i, j) in ALL_TET_PAIRS]


class SemiCornerCuts:
    def __init__(self):
        self.data = [dict(), dict(), dict(), dict()]

    def get(self, i: int, j: int, k: int) -> int:
        return self.data[i].get(_epair(j, k), 0)

    def set(self, i: int, j: int, k: int, v: int):
        self.data[i][_epair(j, k)] = int(v)


class CornerCuts:
    def __init__(self):
        self.data = [0, 0, 0, 0]

    def get(self, i: int) -> int:
        return self.data[i]

    def set(self, i: int, v: int):
        self.data[i] = v

    def __getitem__(self, i: int) -> int:
        return self.data[i]

    def __setitem__(self, i: int, v: int):
        self.data[i] = v


class DiagonalCuts:
    def __init__(self):
        self.data = {}

    def get(self, i: int, j: int) -> int:
        return self.data.get(_epair(i, j), 0)

    def set(self, i: int, j: int, v: int):
        self.data[_epair(i, j)] = v
        for k in range(4):
            if k == i or k == j:
                continue
            for l in range(k + 1, 4):
                if l == i or l == j:
                    continue
                self.data[_epair(k, l)] = v

    def __getitem__(self, ij) -> int:
        return self.get(ij[0], ij[1])


class CombinatorialVertex:
    __slots__ = ("i", "j", "order", "edge_int", "k", "l", "scoop_steiner", "scoop_interior_steiner")

    def __init__(self, i=-1, j=-1, order=-1, edge_int=0, k=-1, l=-1):
        self.i = i
        self.j = j
        self.order = order
        self.edge_int = edge_int
        self.k = k
        self.l = l
        self.scoop_steiner = False
        self.scoop_interior_steiner = False

    def __eq__(self, other) -> bool:
        if self.i == other.i and self.j == other.j and self.order == other.order:
            return True
        if self.i == other.j and self.j == other.i and self.order == self.edge_int - other.order + 1:
            return True
        return False

    def same_edge_as(self, other) -> bool:
        return (self.i == other.i and self.j == other.j) or (self.i == other.j and self.j == other.i)


class DerivedVertex:
    def __init__(self, parents=None, weights=None, dv1=None, dv2=None, t1=None):
        if dv1 is not None:
            # from two derived vertices: v1 * (1 - t1) + v2 * t1
            self.parents = list(dv1.parents) + list(dv2.parents)
            self.weights = [w * (1.0 - t1) for w in dv1.weights] + [w * t1 for w in dv2.weights]
            return
        self.parents = list(parents) if parents is not None else []
        if weights is not None:
            self.weights = list(weights)
        else:
            n = len(self.parents)
            self.weights = [1.0 / n for _ in range(n)] if n else []


class CombFaceType(IntEnum):
    OPEN = 0
    CLOSED_NON_NORMAL = 1
    CLOSED_NORMAL = 2
    TRIANGLE = 3
    QUAD02 = 4
    QUAD13 = 5
    OCTAGON = 6
    SPIRAL = 7
    SMEARED_QUAD = 8
    INTERIOR_HEXAGON = 9
    INTERIOR_CORNER_TYPE = 10
    TUNNEL_QUAD = 11


class CombFace:
    def __init__(self, vertices=None, ftype: CombFaceType = CombFaceType.OPEN):
        self.vertices: List[CombinatorialVertex] = list(vertices) if vertices else []
        self.type = ftype
        self.derived_vertex = DerivedVertex()


class EdgeOccupations:
    def __init__(self, edge_ints: Optional[EdgeInts] = None):
        self.data = {}
        if edge_ints is not None:
            for key, count in edge_ints.data.items():
                self.data[key] = [False] * count

    def _index(self, i, j, index):
        if i > j:
            index = len(self.data[_epair(i, j)]) - 1 - index
        return index

    def get(self, i_or_cv, j=None, index=None) -> bool:
        if j is None:  # CombinatorialVertex
            cv = i_or_cv
            i, j, index = cv.i, cv.j, cv.order - 1
        else:
            i = i_or_cv
        idx = self._index(i, j, index)
        return self.data[_epair(i, j)][idx]

    def set(self, i_or_cv, j=None, index=None, value=None):
        if j is None or (index is None and value is None):
            # set(cv, value)
            cv = i_or_cv
            val = j
            i, j2, index = cv.i, cv.j, cv.order - 1
            idx = self._index(i, j2, index)
            self.data[_epair(i, j2)][idx] = val
        else:
            i = i_or_cv
            idx = self._index(i, j, index)
            self.data[_epair(i, j)][idx] = value

    def to_edge_ints(self, count_occupied: bool) -> EdgeInts:
        result = EdgeInts()
        for key, occ in self.data.items():
            count = 0
            for o in occ:
                if o and count_occupied:
                    count += 1
                if (not o) and (not count_occupied):
                    count += 1
            result.data[key] = count
        return result


# ============================================================================
# Combinatorial signature + TriangleSoup  (from common/triangle_soup.h/.cpp)
# ============================================================================
class CombVertexSigType(IntEnum):
    NORMAL = 0
    SCOOP_FACE_STEINER = 1
    SCOOP_INTERIOR_STEINER = 2
    STANDALONE = 3


@dataclass
class CombVertexKey:
    i: int = -1
    j: int = -1
    order: int = -1
    edge_int: int = 0  # NOT part of identity
    k: int = -1
    type: CombVertexSigType = CombVertexSigType.STANDALONE

    def identity(self) -> Tuple[int, int, int, int, int]:
        return (self.i, self.j, self.order, self.k, int(self.type))


def key_from_cv(cv: CombinatorialVertex) -> CombVertexKey:
    if cv.scoop_steiner:
        typ = CombVertexSigType.SCOOP_FACE_STEINER
        k_val = cv.k
    elif cv.scoop_interior_steiner:
        typ = CombVertexSigType.SCOOP_INTERIOR_STEINER
        k_val = -1
    else:
        typ = CombVertexSigType.NORMAL
        k_val = -1
    i, j, order = cv.i, cv.j, cv.order
    if i > j:
        i, j = j, i
        if typ == CombVertexSigType.NORMAL:
            order = cv.edge_int - order + 1
        else:
            order = cv.edge_int - order
    return CombVertexKey(i, j, order, cv.edge_int, k_val, typ)


def signature_key(s: CombVertexKey):
    """Stable hashable identity for cross-tet merging. Returns None for STANDALONE
    (which must NEVER merge - allocate a fresh vertex instead)."""
    if s.type == CombVertexSigType.STANDALONE:
        return None
    return s.identity()


class TriangleSoup:
    COMB_MERGE = True  # the port always runs true

    def __init__(self):
        self.vertices: List[Vec3] = []
        self.faces: List[List[int]] = []
        self.faces_per_edge: List[List[Tuple[int, int, int]]] = []
        self.dual_positions: List[Vec3] = []
        self.sig_to_idx = {}
        self.signatures: List[CombVertexKey] = []

    def add_local_soup(self, other: "TriangleSoup"):
        if TriangleSoup.COMB_MERGE and other.signatures:
            self.add_local_soup_comb(other)
        else:
            self._add_plain(other.faces, other.vertices)

    def _add_plain(self, local_faces, local_vertices):
        offset = len(self.vertices)
        self.vertices.extend(local_vertices)
        for face in local_faces:
            self.faces.append([idx + offset for idx in face])

    def add_local_soup_comb(self, other: "TriangleSoup"):
        n = len(other.vertices)
        idx_map = [0] * n
        for vi in range(n):
            key = other.signatures[vi]
            if key.type == CombVertexSigType.STANDALONE:
                idx_map[vi] = len(self.vertices)
                self.vertices.append(other.vertices[vi])
                self.signatures.append(key)
            else:
                ident = key.identity()
                if ident in self.sig_to_idx:
                    idx_map[vi] = self.sig_to_idx[ident]
                else:
                    idx_map[vi] = len(self.vertices)
                    self.sig_to_idx[ident] = idx_map[vi]
                    self.vertices.append(other.vertices[vi])
                    self.signatures.append(key)
        for fi in range(len(other.faces)):
            self.faces.append([idx_map[idx] for idx in other.faces[fi]])
            if other.faces_per_edge:
                self.faces_per_edge.append(other.faces_per_edge[fi])
            if other.dual_positions:
                self.dual_positions.append(other.dual_positions[fi])

    def make_vertex_standalone(self, v: int):
        for sig in self.signatures:
            if sig.i == v or sig.j == v:
                self.sig_to_idx.pop(sig.identity(), None)
                sig.type = CombVertexSigType.STANDALONE

    def make_type_standalone(self, target_type: CombVertexSigType):
        for sig in self.signatures:
            if sig.type == target_type:
                self.sig_to_idx.pop(sig.identity(), None)
                sig.type = CombVertexSigType.STANDALONE

    def remap_vertex_indices(self, idx_map: List[int]):
        if not self.signatures:
            return
        for sig in self.signatures:
            if sig.type == CombVertexSigType.STANDALONE:
                continue
            sig.i = idx_map[sig.i]
            sig.j = idx_map[sig.j]
            if sig.i > sig.j:
                sig.i, sig.j = sig.j, sig.i
                if sig.type == CombVertexSigType.NORMAL:
                    sig.order = sig.edge_int - sig.order + 1
                else:
                    sig.order = sig.edge_int - sig.order
            if sig.k >= 0:
                sig.k = idx_map[sig.k]
        for fpe in self.faces_per_edge:
            for t in range(len(fpe)):
                fpe[t] = tuple(idx_map[v] for v in fpe[t])

    def remap_orders(self, edge_isect_ts, occ: EdgeOccupations, count_occupied: bool):
        if not self.signatures:
            return
        for sig in self.signatures:
            if sig.type == CombVertexSigType.STANDALONE:
                continue
            edge_idx = edge_pair_to_index(sig.i, sig.j)
            n = len(edge_isect_ts[edge_idx])
            count = 0
            for k in range(n):
                occ_at_k = occ.get(sig.i, sig.j, k)
                should_count = occ_at_k if count_occupied else (not occ_at_k)
                if should_count:
                    count += 1
                    if count == sig.order:
                        sig.order = k + 1
                        sig.edge_int = n
                        break


# ============================================================================
# NormalCoordinates  (from nc/NC_solver.h/.cpp)
# ============================================================================
class NormalCoordinates:
    def __init__(self, edge_ints: Optional[EdgeInts] = None, compute: bool = True):
        self.edge_ints = EdgeInts(edge_ints) if edge_ints is not None else EdgeInts()
        self.scc = SemiCornerCuts()
        self.cc = CornerCuts()
        self.dc = DiagonalCuts()
        self.scoop = SemiCornerCuts()
        self.triangle_ineq = True
        self.even_sum = True
        if edge_ints is not None and compute:
            self.compute_NC()

    def compute_NC(self):
        self.compute_semi_corner_cuts_and_scoops()
        if (self.even_sum and self.triangle_ineq) or (not self.even_sum):
            scc_reduced = SemiCornerCuts()
            self.compute_corner_cuts_and_reduced_semi_corner_cuts(scc_reduced)
            self.compute_diagonal_cuts(scc_reduced)
        elif self.even_sum:
            reduced_edge_ints = self.scoop_reduction()
            reduced_nc = NormalCoordinates()
            reduced_nc.edge_ints = reduced_edge_ints
            reduced_nc.compute_NC()
            if (not reduced_nc.triangle_ineq) or (not reduced_nc.even_sum):
                raise RuntimeError("compute_NC: scoop reduction did not yield feasible NC!")
            self.cc = reduced_nc.cc
            self.dc = reduced_nc.dc
        else:
            return

    def compute_semi_corner_cuts_and_scoops(self):
        for (i, j, k) in ALL_TET_TRIPLETS:
            feasible, res_ij, res_ik, res_jk = triangle_inequality(
                self.edge_ints[(i, j)], self.edge_ints[(i, k)], self.edge_ints[(j, k)])
            ijk_even_sum = (self.edge_ints[(i, j)] + self.edge_ints[(i, k)] + self.edge_ints[(j, k)]) % 2 == 0
            self.triangle_ineq = feasible and self.triangle_ineq
            self.even_sum = ijk_even_sum and self.even_sum
            if (not ijk_even_sum) and feasible:
                res_ij = 1
                res_ik = 1
                res_jk = 1
            m_ij = self.edge_ints[(i, j)] - res_ij
            m_ik = self.edge_ints[(i, k)] - res_ik
            m_jk = self.edge_ints[(j, k)] - res_jk
            self.scc.set(i, j, k, max(_tdiv2(m_ij + m_ik - m_jk), 0))
            self.scc.set(j, i, k, max(_tdiv2(m_ij + m_jk - m_ik), 0))
            self.scc.set(k, i, j, max(_tdiv2(m_ik + m_jk - m_ij), 0))
            self.scoop.set(i, j, k, _tdiv2(res_jk))
            self.scoop.set(j, i, k, _tdiv2(res_ik))
            self.scoop.set(k, i, j, _tdiv2(res_ij))

    def compute_corner_cuts_and_reduced_semi_corner_cuts(self, scc_reduced: SemiCornerCuts):
        self.cc[0] = min(self.scc.get(0, 1, 2), self.scc.get(0, 1, 3), self.scc.get(0, 2, 3))
        self.cc[1] = min(self.scc.get(1, 0, 2), self.scc.get(1, 0, 3), self.scc.get(1, 2, 3))
        self.cc[2] = min(self.scc.get(2, 0, 1), self.scc.get(2, 0, 3), self.scc.get(2, 1, 3))
        self.cc[3] = min(self.scc.get(3, 0, 1), self.scc.get(3, 0, 2), self.scc.get(3, 1, 2))
        for i in range(4):
            for (j, k) in all_pairs(complement([i])):
                scc_reduced.set(i, j, k, self.scc.get(i, j, k) - self.cc[i])

    def compute_diagonal_cuts(self, scc_reduced: SemiCornerCuts):
        for j in range(1, 4):
            comp = complement([0, j])
            k, l = comp[0], comp[1]
            self.dc.set(0, j, min(scc_reduced.get(k, 0, j), scc_reduced.get(l, 0, j),
                                  scc_reduced.get(j, k, l), scc_reduced.get(0, k, l)))

    def scoop_reduction(self) -> EdgeInts:
        edge_occupations = EdgeOccupations(self.edge_ints)
        for (i, j, k) in ALL_TET_TRIPLETS:
            if self.scoop.get(i, j, k) > 0:
                base, v1, v2 = i, j, k
            elif self.scoop.get(j, k, i) > 0:
                base, v1, v2 = j, k, i
            elif self.scoop.get(k, i, j) > 0:
                base, v1, v2 = k, i, j
            else:
                continue
            scoop_count = self.scoop.get(base, v1, v2)
            for s in range(scoop_count):
                cv_first = CombinatorialVertex(v1, v2, 2 * s + self.scc.get(v1, v2, base) + 1,
                                               self.edge_ints[(v1, v2)], base)
                cv_second = CombinatorialVertex(v1, v2, 2 * s + self.scc.get(v1, v2, base) + 2,
                                                self.edge_ints[(v1, v2)], base)
                self.trace_forever(cv_first, edge_occupations)
                self.trace_forever(cv_second, edge_occupations)
        reduced_edge_ints = EdgeInts()
        for (i, j) in all_pairs([0, 1, 2, 3]):
            occ_count = 0
            for idx in range(self.edge_ints[(i, j)]):
                if not edge_occupations.get(i, j, idx):
                    occ_count += 1
            reduced_edge_ints.set(i, j, occ_count)
        return reduced_edge_ints

    def find_next_vertex(self, v0: CombinatorialVertex) -> CombinatorialVertex:
        if self.even_sum is False:
            raise RuntimeError("combinatorial trace: even_sum is false, cannot trace scoops")
        i, j = v0.i, v0.j
        k = complement([i, j, v0.k])[0]
        order = v0.order
        n_ij = v0.edge_int
        if self.scc.get(i, j, k) + self.scc.get(j, i, k) + 2 * self.scoop.get(k, i, j) != n_ij:
            raise RuntimeError("combinatorial trace: inconsistent scoop and semicorner cut data")
        if order <= self.scc.get(i, j, k):
            return CombinatorialVertex(i, k, order, self.edge_ints.get(i, k), j)
        elif order >= n_ij - self.scc.get(j, i, k) + 1:
            return CombinatorialVertex(j, k, n_ij - order + 1, self.edge_ints.get(j, k), i)
        else:
            if (order - self.scc.get(i, j, k)) % 2 == 0:
                next_order = order - 1
            else:
                next_order = order + 1
            return CombinatorialVertex(i, j, next_order, n_ij, k)

    def trace_forever(self, v0: CombinatorialVertex, edge_occupations: EdgeOccupations):
        cv = v0
        while True:
            if edge_occupations.get(cv):
                return
            edge_occupations.set(cv, True)
            cv = self.find_next_vertex(cv)

    def get_d1_d2(self) -> Tuple[int, int]:
        b1 = 0
        diag_types = 0
        d1 = 0
        d2 = 0
        for i in range(1, 4):
            if self.dc[(b1, i)] > 0:
                diag_types += 1
                if d1 == 0:
                    d1 = self.dc[(b1, i)]
                else:
                    d2 = self.dc[(b1, i)]
        if diag_types > 2:
            raise RuntimeError("get_d1_d2: more than 3 diagonal cut types")
        return d1, d2


# ============================================================================
# boundary curves  (from subgrid_MT/boundary_curve.cpp)
# ============================================================================
def pick_residual_cv(i, j, other, nc: NormalCoordinates, global_inds):
    """Returns (found: bool, cv: CombinatorialVertex)."""
    nij = nc.edge_ints[(i, j)]
    scoop_ij = nc.scoop.get(other, i, j)
    scc_i = nc.scc.get(i, j, other)
    scc_j = nc.scc.get(j, i, other)
    if nij < scoop_ij * 2 + scc_i + scc_j:
        raise RuntimeError("pick_residual_cv: inconsistent scoop and scc counts")
    elif nij == scoop_ij * 2 + scc_i + scc_j:
        return False, CombinatorialVertex()
    elif nij == scoop_ij * 2 + scc_i + scc_j + 1:
        if scoop_ij == 0:
            order = nc.scc.get(i, j, other) + 1
            return True, CombinatorialVertex(i, j, order, nij, other)
        else:
            if scoop_ij % 2 == 0:
                order = nc.scc.get(i, j, other) + scoop_ij + 1
                return True, CombinatorialVertex(i, j, order, nij, other)
            else:
                order = nc.scc.get(i, j, other) + scoop_ij - 1 + 1
                if global_inds[i] > global_inds[j]:
                    order += 2
                return True, CombinatorialVertex(i, j, order, nij, other)
    else:
        raise RuntimeError("pick_residual_cv: more than 1 residual cv on this edge")


def find_next_generic_vertex(cv0: CombinatorialVertex, nc: NormalCoordinates, global_inds):
    i, j = cv0.i, cv0.j
    k = complement([i, j, cv0.k])[0]
    order = cv0.order
    n_ij = cv0.edge_int
    if order <= nc.scc.get(i, j, k):
        return CombinatorialVertex(i, k, order, nc.edge_ints.get(i, k), j)
    elif order >= n_ij - nc.scc.get(j, i, k) + 1:
        return CombinatorialVertex(j, k, n_ij - order + 1, nc.edge_ints.get(j, k), i)
    has_residual, cv_residual = pick_residual_cv(i, j, k, nc, global_inds)
    res_order = cv_residual.order if has_residual else -1
    if order == res_order:
        return cv_residual
    else:
        if order < res_order or (not has_residual):
            if (order - nc.scc.get(i, j, k)) % 2 == 0:
                next_order = order - 1
            else:
                next_order = order + 1
        else:
            if (order - nc.scc.get(i, j, k)) % 2 == 0:
                next_order = order + 1
            else:
                next_order = order - 1
        return CombinatorialVertex(i, j, next_order, n_ij, k)


def trace_generic_curve(start_cv, nc: NormalCoordinates, global_inds, occupations: EdgeOccupations,
                        curve_vertices: List[CombinatorialVertex], keep_occupation_track: bool = True):
    current_cv = start_cv
    if keep_occupation_track:
        occupations.set(current_cv, True)
    curve_vertices.append(current_cv)
    while True:
        next_cv = find_next_generic_vertex(current_cv, nc, global_inds)
        if keep_occupation_track:
            occupations.set(next_cv, True)
        if next_cv == current_cv or next_cv == start_cv:
            break
        curve_vertices.append(next_cv)
        current_cv = next_cv


def pick_scoop_cvs(i, j, other, nc: NormalCoordinates, residual_cv: CombinatorialVertex):
    nij = nc.edge_ints[(i, j)]
    has_residual = residual_cv.i == -1  # default cv (faithful to core's naming quirk)
    start_order = nc.scc.get(i, j, other) + 1
    end_order = nij - nc.scc.get(j, i, other)
    scoop_cvs = []
    for order in range(start_order, end_order + 1):
        cv = CombinatorialVertex(i, j, order, nij, other)
        if (not has_residual) or (not (cv == residual_cv)):
            scoop_cvs.append(cv)
    return scoop_cvs


def pick_residuals_and_scoop_cvs(nc: NormalCoordinates, tet_global_indices):
    residual_cvs = []
    scoop_cvs = []
    for ijk in ALL_TET_TRIPLETS:
        for (i, j) in all_pairs(ijk):
            k = ijk[0] + ijk[1] + ijk[2] - i - j
            found, residual_cv = pick_residual_cv(i, j, k, nc, tet_global_indices)
            if found:
                residual_cvs.append(residual_cv)
            scoop_cvs.extend(pick_scoop_cvs(i, j, k, nc, residual_cv))
    return residual_cvs, scoop_cvs


def boundary_comb_curves(tet_global_indices, edge_isect_counts, populate_normal_curves: bool):
    """Returns (open_curves, scoop_curves, normal_curves, out_edge_occupations)."""
    edge_ints = EdgeInts(list(edge_isect_counts))
    nc = NormalCoordinates()
    nc.edge_ints = edge_ints
    nc.compute_semi_corner_cuts_and_scoops()

    residual_cvs, scoop_cvs = pick_residuals_and_scoop_cvs(nc, tet_global_indices)

    edge_occupations = EdgeOccupations(edge_ints)
    open_curves = []
    for res_cv in residual_cvs:
        if not edge_occupations.get(res_cv):
            curve_vertices = []
            trace_generic_curve(res_cv, nc, tet_global_indices, edge_occupations, curve_vertices)
            open_curves.append(CombFace(curve_vertices, CombFaceType.OPEN))

    scoop_curves = []
    for scoop_cv in scoop_cvs:
        if not edge_occupations.get(scoop_cv):
            curve_vertices = []
            trace_generic_curve(scoop_cv, nc, tet_global_indices, edge_occupations, curve_vertices)
            scoop_curves.append(CombFace(curve_vertices, CombFaceType.CLOSED_NON_NORMAL))

    normal_curves = []
    if populate_normal_curves:
        for (i, j) in ALL_TET_PAIRS:
            for order in range(1, edge_ints[(i, j)] + 1):
                cv = CombinatorialVertex(i, j, order, edge_ints[(i, j)])
                if not edge_occupations.get(cv):
                    curve_vertices = []
                    trace_generic_curve(cv, nc, tet_global_indices, edge_occupations, curve_vertices)
                    normal_curves.append(CombFace(curve_vertices, CombFaceType.CLOSED_NORMAL))
    return open_curves, scoop_curves, normal_curves, edge_occupations


def segment_tet_face(v_a: CombinatorialVertex, v_b: CombinatorialVertex) -> Tuple[int, int, int]:
    if v_a.same_edge_as(v_b):
        return (v_b.i, v_b.j, v_b.k)
    third = v_a.i if (v_a.i != v_b.i and v_a.i != v_b.j) else v_a.j
    return (v_b.i, v_b.j, third)


def curve_segment_faces(comb_face: CombFace) -> List[Tuple[int, int, int]]:
    verts = comb_face.vertices
    is_closed = comb_face.type != CombFaceType.OPEN
    n = len(verts)
    n_segments = n if is_closed else n - 1
    return [segment_tet_face(verts[i], verts[(i + 1) % n]) for i in range(n_segments)]


# ============================================================================
# cv interpolation  (from subgrid_MT/cv_interpolation.cpp)
# ============================================================================
def interpolate_comb_vertex(cv: CombinatorialVertex, tet_positions, edge_isect_ts, bulge: float) -> Vec3:
    i, j, order, n_ij = cv.i, cv.j, cv.order, cv.edge_int
    edge_idx = edge_pair_to_index(i, j)
    if i < j:
        t = edge_isect_ts[edge_idx][order - 1]
    else:
        tmp_order = n_ij - order + 1
        t = 1 - edge_isect_ts[edge_idx][tmp_order - 1]
    pos = (1 - t) * tet_positions[i] + t * tet_positions[j]
    if not cv.scoop_interior_steiner and not cv.scoop_steiner:
        return pos
    if cv.scoop_interior_steiner and cv.scoop_steiner:
        raise ValueError("interpolate_comb_vertex: vertex cannot be both scoop_steiner and scoop_interior_steiner")
    next_order = order + 1
    if i < j:
        next_t = edge_isect_ts[edge_idx][next_order - 1]
    else:
        tmp_next_order = n_ij - next_order + 1
        next_t = 1 - edge_isect_ts[edge_idx][tmp_next_order - 1]
    next_pos = (1 - next_t) * tet_positions[i] + next_t * tet_positions[j]
    mid_pos = 0.5 * (next_pos + pos)
    if cv.scoop_interior_steiner:
        mid_other_edge = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]
                          - tet_positions[i] - tet_positions[j]) / 2.0
        return (1 - bulge) * mid_pos + bulge * mid_other_edge
    else:
        k = cv.k
        face_interior_point = (2.0 * tet_positions[k] + tet_positions[i] + tet_positions[j]) / 4.0
        return (1 - bulge) * mid_pos + bulge * face_interior_point


def interpolate_derived_vertex(dv: DerivedVertex, tet_positions, edge_isect_ts) -> Vec3:
    pos = Vec3.zero()
    for p in range(len(dv.parents)):
        parent_pos = interpolate_comb_vertex(dv.parents[p], tet_positions, edge_isect_ts, -1.0)
        pos = pos + dv.weights[p] * parent_pos
    return pos


def _embed_triangle_face(comb_face, tet_positions, edge_isect_ts) -> TriangleSoup:
    if len(comb_face.vertices) != 3:
        raise ValueError("embed_triangle_face: TRIANGLE must have exactly 3 vertices")
    soup = TriangleSoup()
    for cv in comb_face.vertices:
        soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(cv))
    soup.faces = [[0, 1, 2]]
    return soup


def _embed_quad_face(comb_face, tet_positions, edge_isect_ts) -> TriangleSoup:
    if len(comb_face.vertices) != 4:
        raise ValueError("embed_quad_face: QUAD must have exactly 4 vertices")
    if comb_face.type == CombFaceType.QUAD02:
        tri1, tri2 = [0, 1, 2], [2, 3, 0]
    else:
        tri1, tri2 = [0, 1, 3], [1, 2, 3]
    soup = TriangleSoup()
    for cv in comb_face.vertices:
        soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(cv))
    soup.faces = [tri1, tri2]
    return soup


def _embed_octagon_face(comb_face, tet_positions, edge_isect_ts) -> TriangleSoup:
    if len(comb_face.vertices) != 8:
        raise ValueError("embed_octagon_face: OCTAGON must have exactly 8 vertices")
    n = len(comb_face.vertices)
    soup = TriangleSoup()
    for i in range(n):
        soup.vertices.append(interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, -1.0))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(comb_face.vertices[i]))
        soup.faces.append([n, i, (i + 1) % n])
    soup.vertices.append(interpolate_derived_vertex(comb_face.derived_vertex, tet_positions, edge_isect_ts))
    if TriangleSoup.COMB_MERGE:
        soup.signatures.append(CombVertexKey(-1, -1, -1, 0, -1, CombVertexSigType.STANDALONE))
    return soup


def _embed_spiral_face(comb_face, tet_positions, edge_isect_ts, bulge) -> TriangleSoup:
    n = len(comb_face.vertices)
    soup = TriangleSoup()
    mean_pos = Vec3.zero()
    for i in range(n):
        pos = interpolate_comb_vertex(comb_face.vertices[i], tet_positions, edge_isect_ts, bulge)
        soup.vertices.append(pos)
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(comb_face.vertices[i]))
        mean_pos = mean_pos + pos
        soup.faces.append([n, i, (i + 1) % n])
    mean_pos = mean_pos / float(n)
    soup.vertices.append(mean_pos)
    if TriangleSoup.COMB_MERGE:
        soup.signatures.append(CombVertexKey(-1, -1, -1, 0, -1, CombVertexSigType.STANDALONE))
    return soup


def _embed_interior_hexagon_face(comb_face, tet_positions, edge_isect_ts, bulge) -> TriangleSoup:
    if len(comb_face.vertices) != 9:
        raise ValueError("embed_interior_hexagon_face: INTERIOR_HEXAGON must have exactly 9 vertices")
    soup = TriangleSoup()
    for cv in comb_face.vertices:
        soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(cv))
    soup.faces = [[8, 2, 5], [2, 8, 0], [2, 0, 1], [5, 2, 3], [5, 3, 4], [8, 5, 6], [8, 6, 7]]
    return soup


def _interior_corner_type_bulge_steiner(comb_face, positions, tet_positions, bulge) -> Vec3:
    cv0 = comb_face.vertices[0]
    cv2 = comb_face.vertices[2]
    if cv0.i == cv2.i or cv0.i == cv2.j:
        common_i = cv0.i
    elif cv0.j == cv2.i or cv0.j == cv2.j:
        common_i = cv0.j
    else:
        common_i = -1
    if common_i == -1:
        raise ValueError("first two vertices must be on the same edge")
    l = complement([cv0.i, cv0.j, cv2.i, cv2.j])[0]
    mid_pos = 0.5 * (positions[0] + positions[2])
    mid_other_edge = (tet_positions[l] + tet_positions[common_i]) / 2.0
    return (1 - bulge) * mid_pos + bulge * mid_other_edge


def _embed_interior_corner_type_face(comb_face, tet_positions, edge_isect_ts, mid_scoop_vertices, bulge):
    soup = TriangleSoup()
    for cv in comb_face.vertices:
        soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(cv))
    n = len(soup.vertices)
    if not mid_scoop_vertices:
        for i in range(1, n - 1):
            soup.faces.append([0, i, i + 1])
    else:
        bulge_steiner = _interior_corner_type_bulge_steiner(comb_face, soup.vertices, tet_positions, bulge)
        soup.vertices.append(bulge_steiner)
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(CombVertexKey(-1, -1, -1, 0, -1, CombVertexSigType.STANDALONE))
        for i in range(n):
            soup.faces.append([n, i, (i + 1) % n])
    return soup


def _embed_smeared_face(comb_face, tet_positions, edge_isect_ts, bulge) -> TriangleSoup:
    n = len(comb_face.vertices)
    soup = TriangleSoup()
    for cv in comb_face.vertices:
        soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, bulge))
        if TriangleSoup.COMB_MERGE:
            soup.signatures.append(key_from_cv(cv))
    for i in range(2, n):
        soup.faces.append([1, i, (i + 1) % n])
    return soup


def embed_comb_face(comb_face, tet_positions, edge_isect_ts, mid_scoop_vertices, bulge) -> TriangleSoup:
    if len(comb_face.vertices) < 3:
        raise ValueError("embed_comb_face: CombFace expects at least 3 vertices")
    t = comb_face.type
    if t == CombFaceType.TRIANGLE:
        return _embed_triangle_face(comb_face, tet_positions, edge_isect_ts)
    elif t == CombFaceType.QUAD02 or t == CombFaceType.QUAD13:
        return _embed_quad_face(comb_face, tet_positions, edge_isect_ts)
    elif t == CombFaceType.OCTAGON:
        return _embed_octagon_face(comb_face, tet_positions, edge_isect_ts)
    elif t == CombFaceType.SPIRAL:
        return _embed_spiral_face(comb_face, tet_positions, edge_isect_ts, bulge)
    elif t == CombFaceType.INTERIOR_CORNER_TYPE:
        return _embed_interior_corner_type_face(comb_face, tet_positions, edge_isect_ts, mid_scoop_vertices, bulge)
    elif mid_scoop_vertices and t == CombFaceType.INTERIOR_HEXAGON:
        return _embed_interior_hexagon_face(comb_face, tet_positions, edge_isect_ts, bulge)
    elif t in (CombFaceType.SMEARED_QUAD, CombFaceType.INTERIOR_HEXAGON, CombFaceType.TUNNEL_QUAD):
        return _embed_smeared_face(comb_face, tet_positions, edge_isect_ts, bulge)
    else:
        raise ValueError("embed_comb_face: CombFace type not supported")


# ============================================================================
# simplicial embedding  (from subgrid_MT/simplicial_embedding.cpp)
# ============================================================================
def align_cvs(cv1: CombinatorialVertex, cv2: CombinatorialVertex):
    if not cv1.same_edge_as(cv2):
        raise ValueError("align_cvs: vertices are not on the same edge")
    if cv1.edge_int != cv2.edge_int:
        raise ValueError("align_cvs: vertices do not have the same edge intersection count")
    if cv1.i == cv2.i and cv1.j == cv2.j:
        return cv1, cv2
    elif cv1.i == cv2.j and cv1.j == cv2.i:
        n_ij = cv1.edge_int
        aligned_cv1 = CombinatorialVertex(cv1.i, cv1.j, cv1.order, n_ij, cv1.k)
        aligned_cv2 = CombinatorialVertex(cv1.i, cv1.j, n_ij - cv2.order + 1, n_ij, cv2.k)
        return aligned_cv1, aligned_cv2
    else:
        raise RuntimeError("align_cvs: invalid vertices")


def get_close_scoop_sides(smeared_curve_nc, inner_segment_cv1, inner_segment_cv2) -> List[int]:
    if not inner_segment_cv1.same_edge_as(inner_segment_cv2):
        raise ValueError("get_close_scoop_sides: vertices are not on the same edge")
    i, j = inner_segment_cv1.i, inner_segment_cv1.j
    n_ij = smeared_curve_nc.edge_ints[(i, j)]
    if not (n_ij == inner_segment_cv1.edge_int and n_ij == inner_segment_cv2.edge_int):
        raise ValueError("get_close_scoop_sides: vertices are not on the same edge")
    aligned_cv1, aligned_cv2 = align_cvs(inner_segment_cv1, inner_segment_cv2)
    kl = complement([i, j])
    k, l = kl[0], kl[1]
    smaller_ij_order = min(aligned_cv1.order, aligned_cv2.order)
    close_sides = []
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, k)
            and smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, k)
            and (smaller_ij_order - smeared_curve_nc.scc.get(i, j, k)) % 2 == 1):
        close_sides.append(k)
    if (smaller_ij_order > smeared_curve_nc.scc.get(i, j, l)
            and smaller_ij_order <= n_ij - smeared_curve_nc.scc.get(j, i, l)
            and (smaller_ij_order - smeared_curve_nc.scc.get(i, j, l)) % 2 == 1):
        close_sides.append(l)
    return close_sides


def add_mid_scoop_vertices_spiral(spiral_comb_face: CombFace) -> CombFace:
    new_vertices = []
    num_vertices = len(spiral_comb_face.vertices)
    for idx in range(num_vertices):
        cv = spiral_comb_face.vertices[idx]
        next_cv = spiral_comb_face.vertices[(idx + 1) % num_vertices]
        new_vertices.append(cv)
        if cv.same_edge_as(next_cv):
            if next_cv.k == -1:
                raise RuntimeError("add_mid_scoop_vertices: invalid segment edge with no k index")
            other_vertex = next_cv.k
            aligned_cv, aligned_next_cv = align_cvs(cv, next_cv)
            min_order = min(aligned_cv.order, aligned_next_cv.order)
            mid_cv = CombinatorialVertex(cv.i, cv.j, min_order, cv.edge_int, other_vertex)
            mid_cv.scoop_steiner = True
            new_vertices.append(mid_cv)
    return CombFace(new_vertices, CombFaceType.SPIRAL)


def _add_mid_scoop_vertices_tunnel(tunnel_comb_face, smeared_curve_nc) -> CombFace:
    cv1 = tunnel_comb_face.vertices[0]
    cv2 = tunnel_comb_face.vertices[1]
    close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2)
    if len(close_scoop_sides) < 2:
        return tunnel_comb_face
    elif len(close_scoop_sides) == 2:
        k, l = close_scoop_sides[0], close_scoop_sides[1]
        mid1 = CombinatorialVertex(cv1.i, cv1.j, cv1.order, cv1.edge_int, k)
        mid1.scoop_steiner = True
        mid2 = CombinatorialVertex(cv1.i, cv1.j, cv1.order, cv1.edge_int, l)
        mid2.scoop_steiner = True
        return CombFace([cv1, mid1, cv2, mid2], CombFaceType.TUNNEL_QUAD)
    else:
        raise RuntimeError("add_mid_scoop_vertices: invalid number of close scoop sides")


def _add_mid_scoop_vertices_smeared(even_sum_comb_face, smeared_curve_nc) -> CombFace:
    num_vertices = len(even_sum_comb_face.vertices)
    new_cv_vertices = []
    for idx in range(num_vertices):
        cv1 = even_sum_comb_face.vertices[idx]
        cv2 = even_sum_comb_face.vertices[(idx + 1) % num_vertices]
        new_cv_vertices.append(cv1)
        if not cv1.same_edge_as(cv2):
            continue
        close_scoop_sides = get_close_scoop_sides(smeared_curve_nc, cv1, cv2)
        if len(close_scoop_sides) == 2:
            raise RuntimeError("add_mid_scoop_vertices: invalid scenario should have been handled by TUNNEL_QUAD")
        i, j = cv1.i, cv1.j
        aligned_cv1, aligned_cv2 = align_cvs(cv1, cv2)
        min_order = min(aligned_cv1.order, aligned_cv2.order)
        mid_scoop_cv = CombinatorialVertex(aligned_cv1.i, aligned_cv1.j, min_order, aligned_cv1.edge_int, cv1.k)
        if len(close_scoop_sides) == 0:
            mid_scoop_cv.scoop_interior_steiner = True
            new_cv_vertices.append(mid_scoop_cv)
        elif len(close_scoop_sides) == 1:
            l = close_scoop_sides[0]
            if l == i or l == j:
                raise RuntimeError("add_mid_scoop_vertices: invalid close scoop side")
            mid_scoop_cv.scoop_steiner = True
            mid_scoop_cv.k = l
            new_cv_vertices.append(mid_scoop_cv)
        else:
            raise RuntimeError("add_mid_scoop_vertices: invalid scenario with close scoop sides")
    return CombFace(new_cv_vertices, even_sum_comb_face.type)


def add_mid_scoop_vertices(even_sum_comb_face, smeared_curve_nc) -> CombFace:
    t = even_sum_comb_face.type
    if t == CombFaceType.SPIRAL:
        return add_mid_scoop_vertices_spiral(even_sum_comb_face)
    elif t == CombFaceType.TUNNEL_QUAD:
        return _add_mid_scoop_vertices_tunnel(even_sum_comb_face, smeared_curve_nc)
    elif t in (CombFaceType.SMEARED_QUAD, CombFaceType.INTERIOR_HEXAGON, CombFaceType.INTERIOR_CORNER_TYPE):
        return _add_mid_scoop_vertices_smeared(even_sum_comb_face, smeared_curve_nc)
    else:
        raise RuntimeError("add_mid_scoop_vertices: unsupported comb face type")


# ============================================================================
# evensum construction  (from subgrid_MT/evensum_construction.cpp)
# ============================================================================
def _get_interior_vertex(separated_vs: List[int]) -> int:
    if len(separated_vs) == 4 or len(separated_vs) == 0:
        interior_vertex = -1
    elif len(separated_vs) == 3:
        interior_vertex = complement(separated_vs)[0]
    elif len(separated_vs) == 1:
        interior_vertex = separated_vs[0]
    else:
        raise RuntimeError("get_interior_vertex: Invalid number of separated vertices")
    if interior_vertex < -1 or interior_vertex >= 4:
        raise RuntimeError("get_interior_vertex: computed interior vertex out of range [-1, 3]")
    return interior_vertex


def combface_to_edgeints(comb_face: CombFace) -> EdgeInts:
    result = EdgeInts()
    for cv in comb_face.vertices:
        result.inc(cv.i, cv.j)
    return result


def separated_vertices(cruve_edge_ints: EdgeInts) -> List[int]:
    result = [0]
    for i in range(1, 4):
        if cruve_edge_ints[(0, i)] % 2 == 0:
            result.append(i)
    return result


def _extract_curve_edge_isect_ts(curve_cvs, global_edge_isect_ts):
    curve_edge_ints = EdgeInts(global_edge_isect_ts)
    curve_edge_occupations = EdgeOccupations(curve_edge_ints)
    for cv in curve_cvs:
        curve_edge_occupations.set(cv, True)
    curve_edge_isect_ts = [[] for _ in range(6)]
    for e in range(6):
        i, j = ALL_TET_PAIRS[e]
        for k in range(len(global_edge_isect_ts[e])):
            if curve_edge_occupations.get(i, j, k):
                curve_edge_isect_ts[e].append(global_edge_isect_ts[e][k])
    return curve_edge_isect_ts, curve_edge_occupations


def smeared_quads_construction(smeared_curve_nc, interior_vertex) -> List[CombFace]:
    faces = []
    for ijk in ALL_TET_TRIPLETS:
        for i in ijk:
            j = ijk[1] if ijk[0] == i else ijk[0]
            k = ijk[0] + ijk[1] + ijk[2] - i - j
            n_ij = smeared_curve_nc.edge_ints[(i, j)]
            n_ik = smeared_curve_nc.edge_ints[(i, k)]
            i_is_interior = interior_vertex == i
            for segment_order in range(1, smeared_curve_nc.scc.get(i, j, k) + 1):
                construct = (i_is_interior and segment_order % 2 == 1) or (not i_is_interior and segment_order % 2 == 0)
                if construct and segment_order > 1:
                    cvij_1 = CombinatorialVertex(i, j, segment_order - 1, n_ij, k)
                    cvij_2 = CombinatorialVertex(i, j, segment_order, n_ij, k)
                    cvik_2 = CombinatorialVertex(i, k, segment_order, n_ik, j)
                    cvik_1 = CombinatorialVertex(i, k, segment_order - 1, n_ik, j)
                    faces.append(CombFace([cvij_1, cvij_2, cvik_2, cvik_1], CombFaceType.SMEARED_QUAD))
    return faces


def tunnel_construction(smeared_curve_nc, interior_vertex) -> List[CombFace]:
    faces = []
    for (i, j) in ALL_TET_PAIRS:
        n_ij = smeared_curve_nc.edge_ints[(i, j)]
        for segment_order in range(2, n_ij + 1):
            i_is_interior = interior_vertex == i
            construct = (i_is_interior and segment_order % 2 == 1) or (not i_is_interior and segment_order % 2 == 0)
            if construct:
                cvij_1 = CombinatorialVertex(i, j, segment_order - 1, n_ij)
                cvij_2 = CombinatorialVertex(i, j, segment_order, n_ij)
                faces.append(CombFace([cvij_1, cvij_2], CombFaceType.TUNNEL_QUAD))
    return faces


def interior_region_construction(smeared_curve_nc, interior_vertex) -> List[CombFace]:
    if interior_vertex == -1:
        for ijk in ALL_TET_TRIPLETS:
            i, j, k = ijk[0], ijk[1], ijk[2]
            if smeared_curve_nc.scc.get(i, j, k) % 2 == 1:
                if (smeared_curve_nc.scc.get(i, j, k) != 1 or smeared_curve_nc.scc.get(j, i, k) != 1
                        or smeared_curve_nc.scc.get(k, i, j) != 1):
                    raise RuntimeError("interior_region_construction: invalid contractible curve")
                cv_ij = CombinatorialVertex(i, j, 1, 2, k)
                cv_ik = CombinatorialVertex(i, k, 1, 2, j)
                cv_ji = CombinatorialVertex(j, i, 1, 2, k)
                cv_jk = CombinatorialVertex(j, k, 1, 2, i)
                cv_ki = CombinatorialVertex(k, i, 1, 2, j)
                cv_kj = CombinatorialVertex(k, j, 1, 2, i)
                return [CombFace([cv_ij, cv_ik, cv_ki, cv_kj, cv_jk, cv_ji], CombFaceType.INTERIOR_HEXAGON)]
        return []
    else:
        i = interior_vertex
        jkl = complement([i])
        if smeared_curve_nc.scc.get(i, jkl[0], jkl[1]) == 0:
            j, k = jkl[0], jkl[1]
        elif smeared_curve_nc.scc.get(i, jkl[0], jkl[2]) == 0:
            j, k = jkl[0], jkl[2]
        elif smeared_curve_nc.scc.get(i, jkl[1], jkl[2]) == 0:
            j, k = jkl[1], jkl[2]
        else:
            raise ValueError("interior_region_construction: invalid corner curve")
        l = 6 - i - j - k
        cv_ij = CombinatorialVertex(i, j, 1, smeared_curve_nc.edge_ints[(i, j)], k)
        cv_il = CombinatorialVertex(i, l, 1, smeared_curve_nc.edge_ints[(i, l)], j)
        cv_ik = CombinatorialVertex(i, k, 1, smeared_curve_nc.edge_ints[(i, k)], j)
        cv_kj = CombinatorialVertex(k, j, smeared_curve_nc.scc.get(k, j, i), smeared_curve_nc.edge_ints[(k, j)], i)
        interior_region_vertices = [cv_ij, cv_il, cv_ik, cv_kj]
        n_kj = smeared_curve_nc.edge_ints[(k, j)]
        for scoop_idx in range(smeared_curve_nc.scoop.get(i, j, k)):
            scoop_cv1 = CombinatorialVertex(k, j, 2 * scoop_idx + cv_kj.order + 1, n_kj, i)
            scoop_cv2 = CombinatorialVertex(k, j, 2 * scoop_idx + cv_kj.order + 2, n_kj, i)
            interior_region_vertices.append(scoop_cv1)
            interior_region_vertices.append(scoop_cv2)
        cv_jk = CombinatorialVertex(j, k, smeared_curve_nc.scc.get(j, k, i), smeared_curve_nc.edge_ints[(j, k)], i)
        interior_region_vertices.append(cv_jk)
        return [CombFace(interior_region_vertices, CombFaceType.INTERIOR_CORNER_TYPE)]


def diagonal_surface_construction(tet_positions, global_edge_isect_ts, diagonal_comb_curve,
                                  scoop_mid_vertices, scoop_mid_vertex_bulge) -> TriangleSoup:
    spiral_comb_face = CombFace(list(diagonal_comb_curve.vertices), CombFaceType.SPIRAL)
    if scoop_mid_vertices:
        spiral_comb_face = add_mid_scoop_vertices_spiral(spiral_comb_face)
    return embed_comb_face(spiral_comb_face, tet_positions, global_edge_isect_ts, scoop_mid_vertices, scoop_mid_vertex_bulge)


def nondiagonal_surface_construction(tet_positions, global_edge_isect_ts, nondiagonal_comb_curve,
                                     scoop_mid_vertices, scoop_mid_vertex_bulge) -> TriangleSoup:
    local_edge_isect_ts, curve_edge_occupations = _extract_curve_edge_isect_ts(
        nondiagonal_comb_curve.vertices, global_edge_isect_ts)
    smeared_curve_edge_ints = combface_to_edgeints(nondiagonal_comb_curve)
    interior_vertex = _get_interior_vertex(separated_vertices(smeared_curve_edge_ints))

    smeared_curve_nc = NormalCoordinates(smeared_curve_edge_ints)
    if smeared_curve_nc.triangle_ineq or (not smeared_curve_nc.even_sum):
        raise ValueError("smeared_surface_construction: should satisfy even sum and NOT satisfy triangle inequality")

    comb_faces = []
    comb_faces.extend(smeared_quads_construction(smeared_curve_nc, interior_vertex))
    comb_faces.extend(interior_region_construction(smeared_curve_nc, interior_vertex))
    if scoop_mid_vertices:
        comb_faces.extend(tunnel_construction(smeared_curve_nc, interior_vertex))

    curve_soup = TriangleSoup()
    for face in comb_faces:
        if scoop_mid_vertices:
            face = add_mid_scoop_vertices(face, smeared_curve_nc)
        if len(face.vertices) < 3:
            continue
        curve_soup.add_local_soup(embed_comb_face(face, tet_positions, local_edge_isect_ts,
                                                  scoop_mid_vertices, scoop_mid_vertex_bulge))
    if TriangleSoup.COMB_MERGE:
        curve_soup.remap_orders(global_edge_isect_ts, curve_edge_occupations, True)
    return curve_soup


def evensum_surface_construction(tet_positions, global_edge_isect_ts, non_normal_comb_curves,
                                 scoop_mid_vertices, scoop_mid_vertex_bulge) -> TriangleSoup:
    soup = TriangleSoup()
    for non_normal_comb_curve in non_normal_comb_curves:
        evensum_edge_ints = combface_to_edgeints(non_normal_comb_curve)
        separated_vs = separated_vertices(evensum_edge_ints)
        if len(separated_vs) == 2:
            soup.add_local_soup(diagonal_surface_construction(
                tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge))
        else:
            soup.add_local_soup(nondiagonal_surface_construction(
                tet_positions, global_edge_isect_ts, non_normal_comb_curve, scoop_mid_vertices, scoop_mid_vertex_bulge))
    return soup


# ============================================================================
# normal reconstruction  (from subgrid_MT/normal_reconstruction.cpp)
# ============================================================================
# Construction policy (NOT a runtime option): separate the corner disks from the
# recursive subdivision. This is the shipped behavior and should only be changed
# by a developer who understands the almost-normal recursion. It is deliberately
# a module-local constant rather than a plumbed-through flag.
_NO_CORNER_RECURSION = True


def _shift_orders_past_corners(soup: TriangleSoup, nc: NormalCoordinates):
    for sig in soup.signatures:
        if sig.type == CombVertexSigType.STANDALONE:
            continue
        if sig.i == 4 or sig.j == 4:
            continue
        sig.order += nc.cc[sig.i]
        sig.edge_int = nc.edge_ints[(sig.i, sig.j)]


def corner_construction(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates) -> TriangleSoup:
    soup = TriangleSoup()
    for i in range(4):
        others = complement([i])
        j, k, l = others[0], others[1], others[2]
        for order in range(1, nc.cc[i] + 1):
            triangle_comb_verts = [
                CombinatorialVertex(i, j, order, nc.edge_ints[(i, j)]),
                CombinatorialVertex(i, k, order, nc.edge_ints[(i, k)]),
                CombinatorialVertex(i, l, order, nc.edge_ints[(i, l)]),
            ]
            comb_face = CombFace(triangle_comb_verts, CombFaceType.TRIANGLE)
            soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, False, 0.0))
    return soup


def diagonal_construction(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates) -> TriangleSoup:
    soup = TriangleSoup()
    b1 = 0
    if nc.dc[(b1, 1)] != 0:
        b2 = 1
    elif nc.dc[(b1, 2)] != 0:
        b2 = 2
    elif nc.dc[(b1, 3)] != 0:
        b2 = 3
    else:
        b2 = -1
    if b2 == -1:
        return soup
    head_vs = complement([b1, b2])
    h1, h2 = head_vs[0], head_vs[1]
    diag_count = nc.dc[(b1, b2)]
    b1_cc = nc.cc[b1]
    b2_cc = nc.cc[b2]
    for t in range(diag_count):
        cv0 = CombinatorialVertex(b1, h1, t + b1_cc + 1, nc.edge_ints[(b1, h1)])
        cv1 = CombinatorialVertex(b1, h2, t + b1_cc + 1, nc.edge_ints[(b1, h2)])
        cv2 = CombinatorialVertex(b2, h2, t + b2_cc + 1, nc.edge_ints[(b2, h2)])
        cv3 = CombinatorialVertex(b2, h1, t + b2_cc + 1, nc.edge_ints[(b2, h1)])
        comb_face = CombFace([cv0, cv1, cv2, cv3], CombFaceType.QUAD02)
        soup.add_local_soup(embed_comb_face(comb_face, tet_positions, normal_edge_isect_ts, False, 0.0))
    return soup


def octagon_construction(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates) -> TriangleSoup:
    soup = TriangleSoup()
    b1 = 0
    vd1 = vd2 = vdd = -1
    d1 = d2 = -1
    for other_b1 in (1, 2, 3):
        tmp_d = nc.dc[(b1, other_b1)]
        if tmp_d != 0:
            if vd1 == -1:
                vd1 = other_b1
                d1 = tmp_d
                continue
            else:
                vd2 = other_b1
                d2 = tmp_d
        else:
            vdd = other_b1
    d = d1
    c_b1, c_vd1, c_vd2, c_vdd = nc.cc[b1], nc.cc[vd1], nc.cc[vd2], nc.cc[vdd]
    b1_vdd_mid_vert_1 = CombinatorialVertex(b1, vdd, c_b1 + d, nc.edge_ints[(b1, vdd)])
    b1_vdd_mid_vert_2 = CombinatorialVertex(b1, vdd, c_b1 + d + 1, nc.edge_ints[(b1, vdd)])
    b1_vdd_mid_vertex = DerivedVertex([b1_vdd_mid_vert_1, b1_vdd_mid_vert_2])
    vd1_vd2_mid_vert_1 = CombinatorialVertex(vd1, vd2, c_vd1 + d, nc.edge_ints[(vd1, vd2)])
    vd1_vd2_mid_vert_2 = CombinatorialVertex(vd1, vd2, c_vd1 + d + 1, nc.edge_ints[(vd1, vd2)])
    vd1_vd2_mid_vertex = DerivedVertex([vd1_vd2_mid_vert_1, vd1_vd2_mid_vert_2])
    for t in range(d):
        inner_vertex_t = DerivedVertex(dv1=vd1_vd2_mid_vertex, dv2=b1_vdd_mid_vertex, t1=float(t + 1) / float(d + 1))
        b1vdd_1 = CombinatorialVertex(b1, vdd, t + 1 + c_b1, nc.edge_ints[(b1, vdd)])
        b1vdd_2 = CombinatorialVertex(vdd, b1, t + 1 + c_vdd, nc.edge_ints[(b1, vdd)])
        d1d2_1 = CombinatorialVertex(vd1, vd2, d - t + c_vd1, nc.edge_ints[(vd1, vd2)])
        d1d2_2 = CombinatorialVertex(vd2, vd1, d - t + c_vd2, nc.edge_ints[(vd1, vd2)])
        b1vd1 = CombinatorialVertex(b1, vd1, t + 1 + c_b1, nc.edge_ints[(b1, vd1)])
        b1vd2 = CombinatorialVertex(b1, vd2, t + 1 + c_b1, nc.edge_ints[(b1, vd2)])
        d1vdd = CombinatorialVertex(vd1, vdd, d - t + c_vd1, nc.edge_ints[(vd1, vdd)])
        d2vdd = CombinatorialVertex(vd2, vdd, d - t + c_vd2, nc.edge_ints[(vd2, vdd)])
        octagon_comb_face = CombFace([b1vdd_1, b1vd1, d1d2_1, d1vdd, b1vdd_2, d2vdd, d1d2_2, b1vd2],
                                     CombFaceType.OCTAGON)
        octagon_comb_face.derived_vertex = inner_vertex_t
        soup.add_local_soup(embed_comb_face(octagon_comb_face, tet_positions, normal_edge_isect_ts, False, 0.0))
    return soup


def spiral_construction(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates) -> TriangleSoup:
    i, j = 0, 1
    start_cv = CombinatorialVertex(i, j, nc.cc[i] + 1, nc.edge_ints[(i, j)])
    spiral_comb_verts = []
    dummy_global_inds = [0, 1, 2, 3]
    dummy_edge_occ = EdgeOccupations()
    trace_generic_curve(start_cv, nc, dummy_global_inds, dummy_edge_occ, spiral_comb_verts, False)
    spiral_comb_face = CombFace(spiral_comb_verts, CombFaceType.SPIRAL)
    soup = TriangleSoup()
    soup.add_local_soup(embed_comb_face(spiral_comb_face, tet_positions, normal_edge_isect_ts, False, 0.0))
    return soup


def compute_recursion_center_and_isect_ts(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates):
    def cvp(i, j, base, ei):
        return interpolate_comb_vertex(CombinatorialVertex(i, j, nc.cc[base] + 1, ei), tet_positions,
                                       normal_edge_isect_ts, -1.0)

    p01 = cvp(0, 1, 0, nc.edge_ints[(0, 1)])
    p10 = cvp(1, 0, 1, nc.edge_ints[(0, 1)])
    p02 = cvp(0, 2, 0, nc.edge_ints[(0, 2)])
    p20 = cvp(2, 0, 2, nc.edge_ints[(0, 2)])
    p03 = cvp(0, 3, 0, nc.edge_ints[(0, 3)])
    p30 = cvp(3, 0, 3, nc.edge_ints[(0, 3)])
    p12 = cvp(1, 2, 1, nc.edge_ints[(1, 2)])
    p21 = cvp(2, 1, 2, nc.edge_ints[(1, 2)])
    p13 = cvp(1, 3, 1, nc.edge_ints[(1, 3)])
    p31 = cvp(3, 1, 3, nc.edge_ints[(1, 3)])
    p23 = cvp(2, 3, 2, nc.edge_ints[(2, 3)])
    p32 = cvp(3, 2, 3, nc.edge_ints[(2, 3)])
    subd_center = (p01 + p10 + p02 + p20 + p03 + p30 + p12 + p21 + p13 + p31 + p23 + p32) / 12.0
    hit0, ts0 = segment_triangle_intersection(subd_center, tet_positions[0], p01, p02, p03)
    hit1, ts1 = segment_triangle_intersection(subd_center, tet_positions[1], p10, p12, p13)
    hit2, ts2 = segment_triangle_intersection(subd_center, tet_positions[2], p20, p21, p23)
    hit3, ts3 = segment_triangle_intersection(subd_center, tet_positions[3], p30, p31, p32)
    if not (hit0 and hit1 and hit2 and hit3):
        print("Warning: compute_recursion_center_and_isect_ts: segment-triangle intersection failed; "
              "falling back to greedy recursion intersections")
        ts0 = ts1 = ts2 = ts3 = 1.0
    return subd_center, [ts0, ts1, ts2, ts3]


def recursive_construction(tet_positions, normal_edge_isect_ts, nc: NormalCoordinates,
                           verbose: bool) -> TriangleSoup:
    d1, d2 = nc.get_d1_d2()
    if d1 < d2:
        d1, d2 = d2, d1
    interior_edge_ints = [0, 0, 0, 0]
    b1 = 0
    interior_edge_ints[b1] = d1 - d2 + (0 if _NO_CORNER_RECURSION else nc.cc[b1])
    for other_v in (1, 2, 3):
        if nc.dc[(b1, other_v)] == 0:
            interior_edge_ints[other_v] = 2 * d2 + (0 if _NO_CORNER_RECURSION else nc.cc[other_v])
        elif nc.dc[(b1, other_v)] == d1:
            interior_edge_ints[other_v] = d1 + (0 if _NO_CORNER_RECURSION else nc.cc[other_v])
        elif nc.dc[(b1, other_v)] == d2:
            interior_edge_ints[other_v] = d2 + (0 if _NO_CORNER_RECURSION else nc.cc[other_v])
        else:
            raise RuntimeError("recursive_construction: diagonal cuts cant take any other values")
    new_interior_tet_vertex = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0
    interior_subd_ts_limit = [1.0, 1.0, 1.0, 1.0]
    if _NO_CORNER_RECURSION:
        recursion_center, interior_isect_limits = compute_recursion_center_and_isect_ts(
            tet_positions, normal_edge_isect_ts, nc)
        new_interior_tet_vertex = recursion_center
        for i in range(4):
            interior_subd_ts_limit[i] = interior_isect_limits[i]

    soup = TriangleSoup()
    for triplet in ALL_TET_TRIPLETS:
        i, j, k = triplet[0], triplet[1], triplet[2]
        interior_tet_positions = [new_interior_tet_vertex, tet_positions[i], tet_positions[j], tet_positions[k]]
        interior_tet_edge_isect_ts = [None] * 6
        interior_tet_edge_isect_ts[0] = generate_uniform_isect_ts_single_edge(interior_edge_ints[i], interior_subd_ts_limit[i])
        interior_tet_edge_isect_ts[1] = generate_uniform_isect_ts_single_edge(interior_edge_ints[j], interior_subd_ts_limit[j])
        interior_tet_edge_isect_ts[2] = generate_uniform_isect_ts_single_edge(interior_edge_ints[k], interior_subd_ts_limit[k])
        if _NO_CORNER_RECURSION:
            ts_ij = normal_edge_isect_ts[edge_pair_to_index(i, j)]
            ts_ik = normal_edge_isect_ts[edge_pair_to_index(i, k)]
            ts_jk = normal_edge_isect_ts[edge_pair_to_index(j, k)]
            interior_tet_edge_isect_ts[3] = ts_ij[nc.cc[i]: len(ts_ij) - nc.cc[j]]
            interior_tet_edge_isect_ts[4] = ts_ik[nc.cc[i]: len(ts_ik) - nc.cc[k]]
            interior_tet_edge_isect_ts[5] = ts_jk[nc.cc[j]: len(ts_jk) - nc.cc[k]]
        else:
            interior_tet_edge_isect_ts[3] = normal_edge_isect_ts[edge_pair_to_index(i, j)]
            interior_tet_edge_isect_ts[4] = normal_edge_isect_ts[edge_pair_to_index(i, k)]
            interior_tet_edge_isect_ts[5] = normal_edge_isect_ts[edge_pair_to_index(j, k)]
        sub_soup = normal_surface_construction(interior_tet_positions, interior_tet_edge_isect_ts, verbose)
        if TriangleSoup.COMB_MERGE:
            sub_soup.remap_vertex_indices([4, i, j, k])
            if _NO_CORNER_RECURSION:
                _shift_orders_past_corners(sub_soup, nc)
        soup.add_local_soup(sub_soup)
    if TriangleSoup.COMB_MERGE:
        soup.make_vertex_standalone(4)
    return soup


def normal_surface_construction(tet_positions, normal_edge_isect_ts, verbose: bool) -> TriangleSoup:
    normal_edge_ints = EdgeInts(normal_edge_isect_ts)
    nc = NormalCoordinates(normal_edge_ints)
    if (not nc.triangle_ineq) or (not nc.even_sum):
        raise RuntimeError("Invalid normal coordinates: triangle inequality or even sum condition failed.")
    d1, d2 = nc.get_d1_d2()
    soup = TriangleSoup()
    if d1 == 0 or d2 == 0:
        soup.add_local_soup(corner_construction(tet_positions, normal_edge_isect_ts, nc))
        soup.add_local_soup(diagonal_construction(tet_positions, normal_edge_isect_ts, nc))
    else:
        if d1 == d2 or math.gcd(d1, d2) == 1:
            soup.add_local_soup(corner_construction(tet_positions, normal_edge_isect_ts, nc))
            if d1 == d2:
                soup.add_local_soup(octagon_construction(tet_positions, normal_edge_isect_ts, nc))
            else:
                soup.add_local_soup(spiral_construction(tet_positions, normal_edge_isect_ts, nc))
        else:
            if _NO_CORNER_RECURSION:
                soup.add_local_soup(corner_construction(tet_positions, normal_edge_isect_ts, nc))
            soup.add_local_soup(recursive_construction(tet_positions, normal_edge_isect_ts, nc, verbose))
    return soup


# ============================================================================
# primal entry  (from subgrid_MT/primal_reconstruction.cpp)
# ============================================================================
def unoccupied_edge_isect_ts(edge_isect_ts, edge_occupations: EdgeOccupations):
    result = [[] for _ in range(6)]
    for e_idx in range(6):
        i, j = ALL_TET_PAIRS[e_idx]
        n_ij = len(edge_isect_ts[e_idx])
        for k in range(n_ij):
            if not edge_occupations.get(i, j, k):
                result[e_idx].append(edge_isect_ts[e_idx][k])
    return result


def subgrid_surface(tet_positions, tet_global_indices, edge_isect_ts, open_boundary_curves,
                    scoop_boundary_curves, non_normal_edge_occupations, scoop_mid_vertices,
                    scoop_mid_vertex_bulge, verbose=False) -> TriangleSoup:
    soup = TriangleSoup()
    # open curves are discarded
    evensum_soup = evensum_surface_construction(tet_positions, edge_isect_ts, scoop_boundary_curves,
                                                scoop_mid_vertices, scoop_mid_vertex_bulge)
    soup.add_local_soup(evensum_soup)
    normal_edge_isect_ts = unoccupied_edge_isect_ts(edge_isect_ts, non_normal_edge_occupations)
    normal_soup = normal_surface_construction(tet_positions, normal_edge_isect_ts, verbose)
    if TriangleSoup.COMB_MERGE:
        normal_soup.remap_orders(edge_isect_ts, non_normal_edge_occupations, False)
    soup.add_local_soup(normal_soup)
    if TriangleSoup.COMB_MERGE:
        idx_map = [int(tet_global_indices[0]), int(tet_global_indices[1]),
                   int(tet_global_indices[2]), int(tet_global_indices[3])]
        soup.remap_vertex_indices(idx_map)
        soup.make_type_standalone(CombVertexSigType.SCOOP_INTERIOR_STEINER)
    return soup


# ============================================================================
# Public output struct + folded entry point
# ============================================================================
@dataclass
class PrimalResult:
    vertices: List[Vec3] = field(default_factory=list)      # list of Vec3 (also iterable as (x,y,z))
    faces: List[List[int]] = field(default_factory=list)    # polygons (n-gons), indices into vertices
    signatures: List[CombVertexKey] = field(default_factory=list)  # parallel to vertices; use signature_key()
    non_even: bool = False                                   # True if this tet is NOT even-sum (open curves present)


def subgrid_primal(tet_positions, tet_global_indices, edge_isect_ts,
                   scoop_mid_vertices: bool = True, scoop_mid_vertex_bulge: float = 1e-4) -> PrimalResult:
    """Folds boundary_comb_curves inside; input matches construction_policy.md "Input".

    tet_positions:       list of 4 points (Vec3 or (x,y,z))
    tet_global_indices:  list of 4 ints
    edge_isect_ts:       list of 6 lists of floats, edges (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    """
    TriangleSoup.COMB_MERGE = True
    pos = [_as_vec3(p) for p in tet_positions]
    idx = list(tet_global_indices)
    counts = [len(edge_isect_ts[e]) for e in range(6)]

    open_curves, scoop_curves, _normal_curves, non_normal_edge_occupations = boundary_comb_curves(
        idx, counts, populate_normal_curves=False)

    soup = subgrid_surface(pos, idx, edge_isect_ts, open_curves, scoop_curves,
                           non_normal_edge_occupations, scoop_mid_vertices, scoop_mid_vertex_bulge, False)

    return PrimalResult(vertices=soup.vertices, faces=soup.faces, signatures=soup.signatures,
                        non_even=bool(open_curves))


# ============================================================================
# Dual: QEF solver  (from subgrid_MT/dual/qef_solver.cpp; Eigen ldlt -> hand-rolled 3x3)
# ============================================================================
def _solve3x3(A, b: Vec3) -> Vec3:
    """Solve A x = b for a 3x3 matrix via analytic inverse. A is SPD here (sum of
    n n^T + reg*I), so this stands in for Eigen's A.ldlt().solve(b)."""
    c00 = A[1][1] * A[2][2] - A[1][2] * A[2][1]
    c01 = A[1][0] * A[2][2] - A[1][2] * A[2][0]
    c02 = A[1][0] * A[2][1] - A[1][1] * A[2][0]
    det = A[0][0] * c00 - A[0][1] * c01 + A[0][2] * c02
    i00 = c00 / det
    i01 = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det
    i02 = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det
    i10 = -c01 / det
    i11 = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det
    i12 = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det
    i20 = c02 / det
    i21 = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det
    i22 = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det
    return Vec3(i00 * b.x + i01 * b.y + i02 * b.z, i10 * b.x + i11 * b.y + i12 * b.z,
                i20 * b.x + i21 * b.y + i22 * b.z)


def QEF_from_boundary_polygons(polygons, positions, normals, regularizer_weight=1e-3):
    dual_points = []
    for poly in polygons:
        if len(poly) < 2:
            dual_points.append(Vec3.undefined())
            continue
        avg_pos = Vec3.zero()
        for vidx in poly:
            avg_pos = avg_pos + positions[vidx]
        avg_pos = avg_pos / float(len(poly))
        A = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        b = Vec3.zero()
        for vidx in poly:
            n = normals[vidx]
            p = positions[vidx]
            A[0][0] += n.x * n.x; A[0][1] += n.x * n.y; A[0][2] += n.x * n.z
            A[1][0] += n.y * n.x; A[1][1] += n.y * n.y; A[1][2] += n.y * n.z
            A[2][0] += n.z * n.x; A[2][1] += n.z * n.y; A[2][2] += n.z * n.z
            b = b + n * dot(n, p)
        A[0][0] += regularizer_weight
        A[1][1] += regularizer_weight
        A[2][2] += regularizer_weight
        b = b + regularizer_weight * avg_pos
        dual_points.append(_solve3x3(A, b))
    return dual_points


def get_centroids(polygons, positions):
    centroids = []
    for poly in polygons:
        centroid = Vec3.zero()
        for vidx in poly:
            centroid = centroid + positions[vidx]
        centroid = centroid / float(len(poly))
        centroids.append(centroid)
    return centroids


def projected_dual_position(dual_p: Vec3, centroid: Vec3, tet_positions) -> Vec3:
    if is_inside_tet(dual_p, tet_positions):
        return dual_p
    eps = 1e-6
    tet_center = (tet_positions[0] + tet_positions[1] + tet_positions[2] + tet_positions[3]) / 4.0
    adjusted_centroid = centroid + 1e-6 * (tet_center - centroid)
    _hit, isect_t = segment_tet_intersection(dual_p, adjusted_centroid, tet_positions)
    dir = adjusted_centroid - dual_p
    return dual_p + (isect_t + eps) * dir


def projected_dual_positions(dual_positions, centroids, tet_positions):
    out = []
    for idx in range(len(dual_positions)):
        dp = dual_positions[idx]
        if not dp.is_defined():
            out.append(dp)
            continue
        out.append(projected_dual_position(dp, centroids[idx], tet_positions))
    return out


def dual_positions_from_isect_data(tet_positions, boundary_polygons, isect_positions, isect_normals,
                                   regularizer_weight, projected):
    dual_positions = QEF_from_boundary_polygons(boundary_polygons, isect_positions, isect_normals, regularizer_weight)
    if projected:
        centroids = get_centroids(boundary_polygons, isect_positions)
        dual_positions = projected_dual_positions(dual_positions, centroids, tet_positions)
    return dual_positions


# ============================================================================
# Dual construction  (from subgrid_MT/dual/dual_construction.cpp)
# ============================================================================
def normals_from_cv_polygons(curves, edge_isect_ts, edge_isect_normals):
    normals = []
    for curve in curves:
        for cv in curve.vertices:
            i, j = cv.i, cv.j
            edge_idx = edge_pair_to_index(i, j)
            idx = (cv.order - 1) if i < j else (cv.edge_int - cv.order)
            normals.append(_as_vec3(edge_isect_normals[edge_idx][idx]))
    return normals


def build_dual_local_soup(curves, tet_positions, edge_isect_ts) -> TriangleSoup:
    soup = TriangleSoup()
    vert_offset = 0
    for curve in curves:
        n = len(curve.vertices)
        polygon = [0] * n
        fpe = curve_segment_faces(curve)
        for vi in range(n):
            cv = curve.vertices[vi]
            soup.vertices.append(interpolate_comb_vertex(cv, tet_positions, edge_isect_ts, -1.0))
            soup.signatures.append(key_from_cv(cv))
            polygon[vi] = vert_offset + vi
        soup.faces.append(polygon)
        soup.faces_per_edge.append(fpe)
        vert_offset += n
    return soup


def dual_subgrid_surface(tet_positions, tet_indices, edge_isect_ts, edge_isect_normals,
                         open_curves, scoop_curves, normal_curves, reg_alpha, project_duals) -> TriangleSoup:
    curves = list(scoop_curves) + list(normal_curves)
    if not curves:
        return TriangleSoup()
    normals = normals_from_cv_polygons(curves, edge_isect_ts, edge_isect_normals)
    soup = build_dual_local_soup(curves, tet_positions, edge_isect_ts)
    soup.dual_positions = dual_positions_from_isect_data(
        tet_positions, soup.faces, soup.vertices, normals, reg_alpha, project_duals)
    idx_map = [int(tet_indices[0]), int(tet_indices[1]), int(tet_indices[2]), int(tet_indices[3])]
    soup.remap_vertex_indices(idx_map)
    return soup


# ============================================================================
# Dual entry (folded)
# ============================================================================
@dataclass
class DualResult:
    vertices: List[Vec3] = field(default_factory=list)
    faces: List[List[int]] = field(default_factory=list)                 # boundary loops
    faces_per_edge: List[List[Tuple[int, int, int]]] = field(default_factory=list)  # global {i,j,k} per polygon edge
    dual_positions: List[Vec3] = field(default_factory=list)             # one QEF dual point per polygon
    signatures: List[CombVertexKey] = field(default_factory=list)        # all NORMAL type
    non_even: bool = False
    non_normal: bool = False


def subgrid_dual(tet_positions, tet_global_indices, edge_isect_ts, edge_isect_normals,
                 reg_alpha: float = 0.1, project_duals: bool = False) -> DualResult:
    """Folds boundary_comb_curves inside; input matches construction_policy.md "Input".

    edge_isect_normals: list of 6 lists of surface normals (Vec3 or (x,y,z)),
    parallel to edge_isect_ts.
    """
    TriangleSoup.COMB_MERGE = True
    pos = [_as_vec3(p) for p in tet_positions]
    idx = list(tet_global_indices)
    counts = [len(edge_isect_ts[e]) for e in range(6)]

    open_curves, scoop_curves, normal_curves, _occ = boundary_comb_curves(
        idx, counts, populate_normal_curves=True)

    soup = dual_subgrid_surface(pos, idx, edge_isect_ts, edge_isect_normals,
                                open_curves, scoop_curves, normal_curves, reg_alpha, project_duals)

    return DualResult(vertices=soup.vertices, faces=soup.faces, faces_per_edge=soup.faces_per_edge,
                      dual_positions=soup.dual_positions, signatures=soup.signatures,
                      non_even=bool(open_curves),
                      non_normal=not (len(open_curves) == 0 and len(scoop_curves) == 0))
