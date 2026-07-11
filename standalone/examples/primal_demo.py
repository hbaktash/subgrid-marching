#!/usr/bin/env python3
"""Minimal demo + global-merge snippet for the standalone primal port.
    python3 primal_demo.py
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import subgrid_mt as smt  # noqa: E402


def main():
    # One tet with a few intersections on its edges.
    pos = [smt.Vec3(0, 0, 0), smt.Vec3(1, 0, 0), smt.Vec3(0, 1, 0), smt.Vec3(0, 0, 1)]
    gidx = [0, 1, 2, 3]
    ts = [[] for _ in range(6)]  # edges (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    # A single loop cutting off vertex 0 (crosses the 3 edges incident to it):
    # every face is even-sum, so this is a normal curve -> one triangle.
    ts[0] = [0.3]
    ts[1] = [0.3]
    ts[2] = [0.3]

    r = smt.subgrid_primal(pos, gidx, ts)
    print(f"vertices: {len(r.vertices)}, faces: {len(r.faces)}, non_even: {r.non_even}")

    # ---- Global cross-tet merge sketch (mirrors add_local_soup_comb) ----
    # Accumulate a PrimalResult for every active tet, then merge by signature.
    global_v = []
    global_f = []
    key_to_gid = {}

    def add_result(res):
        local_to_global = [0] * len(res.vertices)
        for vi in range(len(res.vertices)):
            key = smt.signature_key(res.signatures[vi])
            if key is None:  # STANDALONE -> always fresh
                local_to_global[vi] = len(global_v)
                global_v.append(res.vertices[vi])
            elif key in key_to_gid:
                local_to_global[vi] = key_to_gid[key]
            else:
                local_to_global[vi] = len(global_v)
                key_to_gid[key] = len(global_v)
                global_v.append(res.vertices[vi])
        for f in res.faces:
            global_f.append([local_to_global[idx] for idx in f])

    add_result(r)
    print(f"merged global vertices: {len(global_v)}, faces: {len(global_f)}")


if __name__ == "__main__":
    main()
