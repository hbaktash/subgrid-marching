#!/usr/bin/env python3
"""Minimal demo + global-merge snippet for the standalone dual port.
    python3 dual_demo.py
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import subgrid_mt as smt  # noqa: E402


def main():
    # One tet with a few intersections, plus a surface normal per intersection
    # (the dual QEF solve needs them; here we just make some up).
    pos = [smt.Vec3(0, 0, 0), smt.Vec3(1, 0, 0), smt.Vec3(0, 1, 0), smt.Vec3(0, 0, 1)]
    gidx = [0, 1, 2, 3]
    ts = [[] for _ in range(6)]    # edges (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    nrm = [[] for _ in range(6)]   # parallel to ts
    # A single loop cutting off vertex 0 (crosses the 3 edges incident to it).
    ts[0] = [0.3]; nrm[0] = [(1, 1, 1)]
    ts[1] = [0.3]; nrm[1] = [(1, 1, 1)]
    ts[2] = [0.3]; nrm[2] = [(1, 1, 1)]

    r = smt.subgrid_dual(pos, gidx, ts, nrm, reg_alpha=1e-3, project_duals=False)
    print(f"boundary-loop vertices: {len(r.vertices)}, loops (faces): {len(r.faces)}, "
          f"dual points: {len(r.dual_positions)}, non_normal: {r.non_normal}")

    # ---- Global cross-tet merge sketch (mirrors add_local_soup) ----
    # Boundary-loop vertices merge across tets exactly like primal (signature_key).
    # r.dual_positions[f] is the QEF point for loop f and r.faces_per_edge[f] gives
    # the tet face {i,j,k} each loop edge crosses; the downstream dual-contouring
    # step uses those two to stitch dual points into the final mesh.
    global_v = []
    key_to_gid = {}
    local_to_global = [0] * len(r.vertices)
    for vi in range(len(r.vertices)):
        key = smt.signature_key(r.signatures[vi])
        if key is None:  # STANDALONE -> always fresh
            local_to_global[vi] = len(global_v)
            global_v.append(r.vertices[vi])
        elif key in key_to_gid:
            local_to_global[vi] = key_to_gid[key]
        else:
            local_to_global[vi] = len(global_v)
            key_to_gid[key] = len(global_v)
            global_v.append(r.vertices[vi])
    print(f"merged global boundary vertices: {len(global_v)}")


if __name__ == "__main__":
    main()
