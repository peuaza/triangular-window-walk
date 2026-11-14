"""The k-sweep table for the merged paper: meeting fraction by k, and where the gain stops.

Reads probe_colour_k/k_cells.json and writes _k_table.tex.
"""
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "k_cells.json"), encoding="utf-8"))
C = D["cells"]
RULES = ["pattern", "runs", "peak"]
NS = ["512", "1024", "2048"]


def nm(r):
    return r"\texttt{%s}" % r


def f(x, nd=3):
    return "--" if x is None else ("%%.%df" % nd) % x


def main():
    L = [r"\begin{tabular}{@{}llrrrrrrr@{}}", r"\toprule",
         r" & & \multicolumn{4}{c}{meeting fraction} & \multicolumn{3}{c}{step gain (tapes of 300)} \\",
         r"\cmidrule(lr){3-6}\cmidrule(lr){7-9}",
         r"rule & $n$ & $k=1$ & $k=2$ & $k=3$ & $k=4$ & $1\to2$ & $2\to3$ & $3\to4$ \\",
         r"\midrule"]
    for r in RULES:
        for n in NS:
            vals = []
            for k in ("1", "2", "3", "4"):
                c = C[r][k].get(n)
                vals.append(c["meet"] if c else None)
            if any(v is None for v in vals):
                continue
            fr = [v / 300.0 for v in vals]
            d = [vals[1] - vals[0], vals[2] - vals[1], vals[3] - vals[2]]
            L.append(" & ".join([nm(r) if n == "512" else "", n,
                                 f(fr[0]), f(fr[1]), f(fr[2]), f(fr[3]),
                                 "%+d" % d[0], "%+d" % d[1], "%+d" % d[2]]) + r" \\")
        L.append(r"\midrule")
    L = L[:-1]
    L += [r"\bottomrule", r"\end{tabular}"]
    out = os.path.join(HERE, "_k_table.tex")
    io.open(out, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
