#!/usr/bin/env python3
"""Verify the Python primal port against fixtures dumped from the real core by
the C++ `verify_port` harness (ground truth = src/subgrid_MT core).

Usage:
    python3 verify_primal.py [fixtures_primal.txt]
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import subgrid_mt as smt  # noqa: E402

POS_TOL = 1e-9


def read_fixtures(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    i = 0
    n = len(lines)
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
        non_even = bool(int(lines[i].split()[1]))
        i += 1
        nv = int(lines[i].split()[1])
        i += 1
        V = []
        for _ in range(nv):
            x, y, z = lines[i].split()
            V.append((float(x), float(y), float(z)))
            i += 1
        nf = int(lines[i].split()[1])
        i += 1
        F = []
        for _ in range(nf):
            parts = [int(x) for x in lines[i].split()]
            F.append(parts[1:])
            i += 1
        ns = int(lines[i].split()[1])
        i += 1
        S = []
        for _ in range(ns):
            a = [int(x) for x in lines[i].split()]
            S.append(tuple(a))  # (i, j, order, edge_int, k, type)
            i += 1
        cases.append((counts, idx, ts, non_even, V, F, S))
    return cases


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "fixtures_primal.txt")
    cases = read_fixtures(path)

    compared = 0
    mismatches = 0
    max_pos_diff = 0.0
    pos = [smt.Vec3(0, 0, 0), smt.Vec3(1, 0, 0), smt.Vec3(0, 1, 0), smt.Vec3(0, 0, 1)]

    for (counts, idx, ts, non_even, V, F, S) in cases:
        try:
            r = smt.subgrid_primal(pos, idx, ts)
        except Exception as e:  # noqa: BLE001
            mismatches += 1
            if mismatches <= 20:
                print(f"PORT THREW on edges={counts}: {e}")
            continue
        compared += 1

        ok = True
        if r.non_even != non_even:
            ok = False
        if len(r.vertices) != len(V) or len(r.faces) != len(F) or len(r.signatures) != len(S):
            ok = False
        if ok:
            for a, b in zip(r.vertices, V):
                max_pos_diff = max(max_pos_diff, abs(a.x - b[0]), abs(a.y - b[1]), abs(a.z - b[2]))
            for a, b in zip(r.faces, F):
                if list(a) != list(b):
                    ok = False
                    break
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
                      f"portF={len(r.faces)} refF={len(F)} ne(p/r)={r.non_even}/{non_even}")

    print("\n=== python primal verification ===")
    print(f"cases:            {len(cases)}")
    print(f"compared:         {compared}")
    print(f"mismatches:       {mismatches}")
    print(f"max position diff:{max_pos_diff:.3e}")
    print("PASS" if mismatches == 0 else "FAIL")
    return 0 if mismatches == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
