"""Attack the double-lattice rule before believing anything about it.

phi(w) = P_L(w) xor P_R(w),  P_L = parity of w[0], w[D], w[2D], ...,
                            P_R = parity of w[m-1], w[m-1-D], ...
D = D_of(m) = max(2, round(log2 m)).

The walk is lf_lib's action table, unmodified.  Three charge readings are computed from the same
trace, so the comparison is within one walk:
  A  read both lattices per query           2*ceil(m/D)
  B  read only the lattice that a move shifts  ceil(m/D)   (the attack)
  C  full window scan per query             m
Reported per n: queries, and e, g, G for each reading, plus the spread g/e.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, "probe_log_family"))
import lf_lib as LF  # noqa: E402
sys.path.insert(0, os.path.join(PAPER, "probe_triple_nat"))
import tr_lib as T  # noqa: E402


def make_phi_dl():
    def phi(sg, ell, m):
        D = LF.D_of(m)
        pl = int(np.asarray(sg[ell:ell + m:D]).sum() & 1)
        pr = int(np.asarray(sg[ell + m - 1:ell - 1 if ell else None:-D]).sum() & 1)
        return pl ^ pr
    return phi


def costs_from_trace(tr):
    """tr has 'm' per query (list).  Return the three charges."""
    ms = np.asarray(tr["m"], np.int64)
    Ds = np.array([LF.D_of(int(m)) for m in ms], np.int64)
    lat = -(-ms // Ds)                      # ceil(m/D)
    return {"A": int(2 * lat.sum()), "B": int(lat.sum()), "C": int(ms.sum())}


def main(ns):
    phi = make_phi_dl()
    out = {}
    for n in ns:
        rows = {"A": [], "B": [], "C": []}
        qs = []
        for i, a, s, X in T.tapes(n):
            r = LF.walk(a, n, phi, ks=LF.KSHIFT)
            if r.get("end") != "meet":
                continue
            tr = r["tr"]
            c = costs_from_trace(tr)
            qs.append(len(tr["m"]))
            for k in rows:
                rows[k].append(c[k])
        res = {"queries_mean": float(np.mean(qs)), "tapes": len(qs)}
        for k, v in rows.items():
            v = np.asarray(v, float)
            ge = np.exp(np.mean(np.log(np.maximum(v, 1.0))))
            res[k] = {"e": float(v.mean()), "g": float(v.max()), "G": float(ge),
                      "g_over_e": float(v.max() / v.mean()), "G_over_e": float(ge / v.mean()),
                      "e_sd": float(v.std())}
        out[n] = res
        print(n, "tapes", res["tapes"], "q", round(res["queries_mean"], 1), flush=True)
        for k in ("A", "B", "C"):
            d = res[k]
            print("   %s  e=%9.1f  g/e=%.4f  G/e=%.4f  sd=%8.1f" % (k, d["e"], d["g_over_e"], d["G_over_e"], d["e_sd"]), flush=True)
    import json
    json.dump(out, open(os.path.join(HERE, "dl_attack.json"), "w"), indent=1)
    print("wrote dl_attack.json")


if __name__ == "__main__":
    main([int(x) for x in (sys.argv[1:] or ["256", "512", "1024"])])
