#!/usr/bin/env python3
"""Verify the Python dual port against fixtures dumped from the real core by the
C++ `verify_port_dual` harness (ground truth = src/subgrid_MT/dual core).

Usage:
    python3 verify_dual.py [fixtures_dual.txt] [project]
"""
import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import subgrid_mt as smt  # noqa: E402

POS_TOL = 1e-9
DUAL_TOL = 1e-6  # QEF solve: analytic inverse vs Eigen LDLT


def read_fixtures(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    i, n = 0, len(lines)
    cases = []
    while i < n:
        if not lines[i].startswith("C "):
            i += 1
            continue
        counts = [int(x) for x in lines[i].split()[1:]]
        i += 1
        idx = [int(x) for x in lines[i].split()[1:]]
        i += 1
        ts = []
        for _ in range(6):
            parts = lines[i].split()
            cnt = int(parts[1])
            ts.append([float(x) for x in parts[2:2 + cnt]])
            i += 1
        nrm = []
        for _ in range(6):
            parts = lines[i].split()
            cnt = int(parts[1])
            vals = [float(x) for x in parts[2:2 + 3 * cnt]]
            nrm.append([(vals[3 * k], vals[3 * k + 1], vals[3 * k + 2]) for k in range(cnt)])
            i += 1
        nv = int(lines[i].split()[1]); i += 1
        V = []
        for _ in range(nv):
            x, y, z = lines[i].split(); V.append((float(x), float(y), float(z))); i += 1
        nf = int(lines[i].split()[1]); i += 1
        F = []
        for _ in range(nf):
            parts = [int(x) for x in lines[i].split()]; F.append(parts[1:]); i += 1
        ne = int(lines[i].split()[1]); i += 1
        E = []
        for _ in range(ne):
            parts = [int(x) for x in lines[i].split()]
            cnt = parts[0]
            tris = [(parts[1 + 3 * k], parts[2 + 3 * k], parts[3 + 3 * k]) for k in range(cnt)]
            E.append(tris); i += 1
        nd = int(lines[i].split()[1]); i += 1
        D = []
        for _ in range(nd):
            parts = lines[i].split()
            if parts[0] == "nan":
                D.append(None)
            else:
                D.append((float(parts[0]), float(parts[1]), float(parts[2])))
            i += 1
        ns = int(lines[i].split()[1]); i += 1
        S = []
        for _ in range(ns):
            a = [int(x) for x in lines[i].split()]
            S.append(tuple(a))  # (i, j, order, edge_int, k, type)
            i += 1
        cases.append((counts, idx, ts, nrm, V, F, E, D, S))
    return cases


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "fixtures_dual.txt")
    project = len(sys.argv) > 2 and sys.argv[2] == "project"
    cases = read_fixtures(path)

    compared = 0
    mismatches = 0
    max_pos_diff = 0.0
    max_dual_diff = 0.0
    pos = [smt.Vec3(0, 0, 0), smt.Vec3(1, 0, 0), smt.Vec3(0, 1, 0), smt.Vec3(0, 0, 1)]

    for (counts, idx, ts, nrm, V, F, E, D, S) in cases:
        try:
            r = smt.subgrid_dual(pos, idx, ts, nrm, 0.1, project)
        except Exception as e:  # noqa: BLE001
            mismatches += 1
            if mismatches <= 20:
                print(f"PORT THREW on edges={counts}: {e}")
            continue
        compared += 1

        ok = True
        if (len(r.vertices) != len(V) or len(r.faces) != len(F) or len(r.faces_per_edge) != len(E)
                or len(r.dual_positions) != len(D) or len(r.signatures) != len(S)):
            ok = False
        if ok:
            for a, b in zip(r.vertices, V):
                max_pos_diff = max(max_pos_diff, abs(a.x - b[0]), abs(a.y - b[1]), abs(a.z - b[2]))
            for a, b in zip(r.faces, F):
                if list(a) != list(b):
                    ok = False
                    break
        if ok:
            for a, b in zip(r.faces_per_edge, E):
                if [tuple(t) for t in a] != [tuple(t) for t in b]:
                    ok = False
                    break
        if ok:
            for a, b in zip(r.dual_positions, D):
                ad = a.is_defined()
                bd = b is not None
                if ad != bd:
                    ok = False
                    break
                if ad:
                    max_dual_diff = max(max_dual_diff, abs(a.x - b[0]), abs(a.y - b[1]), abs(a.z - b[2]))
        if ok:
            for sig, b in zip(r.signatures, S):
                got = (sig.i, sig.j, sig.order, sig.edge_int, sig.k, int(sig.type))
                if got != b:
                    ok = False
                    break
        if not ok:
            mismatches += 1
            if mismatches <= 20:
                print(f"MISMATCH edges={counts}: portV={len(r.vertices)} refV={len(V)} "
                      f"portF={len(r.faces)} refF={len(F)}")

    print("\n=== python dual verification (project=%d) ===" % int(project))
    print(f"cases:             {len(cases)}")
    print(f"compared:          {compared}")
    print(f"mismatches:        {mismatches}")
    print(f"max position diff: {max_pos_diff:.3e}")
    print(f"max dual pos diff: {max_dual_diff:.3e}")
    print("PASS" if mismatches == 0 else "FAIL")
    return 0 if mismatches == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
