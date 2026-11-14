"""Two more tables for the BTHM-C paper, generated from probe_cycle/cyc_cells.json.

_rep_table.tex : single trim versus double trim on the repeat, per rule at n = 512.
_agg_table.tex : the aggregate at the three sizes (over the twenty rules that were run).
"""
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "cyc_cells.json"), encoding="utf-8"))
R = D["rules"]
TS = D["trim_sensitivity"]
ORDER = ["canonical", "runs", "dense2", "peak", "trip", "oddrun", "pattern", "LAS", "MAXP", "KCHG",
         "BAL>m/2", "BAL>m/4", "NBLK=1", "NBLK>2", "NBLK>3", "AMAX", "AMIN", "longrun", "pal", "hashmod"]


def f(x, nd=3):
    return "--" if x is None else ("%%.%df" % nd) % x


def i0(x):
    return "--" if x is None else "%d" % round(x)


def tname(r):
    return r"\texttt{%s}" % r.replace(">", r"$>$").replace("=", r"$=$").replace("_", r"\_")


def rep_table():
    L = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{meet} & \multicolumn{2}{c}{reads (median)} & \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
         r"rule & left & right & left & right & $\Delta$ \\", r"\midrule"]
    for r in ORDER:
        v = TS["512"].get(r)
        if not v:
            continue
        d = (v["trim2_meet_frac"] - v["trim1_meet_frac"])
        L.append(" & ".join([tname(r), f(v["trim1_meet_frac"]), f(v["trim2_meet_frac"]),
                             i0(v["trim1_reads_med_all"]), i0(v["trim2_reads_med_all"]),
                             ("%+.3f" % d)]) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    io.open(os.path.join(HERE, "_rep_table.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


def agg_table():
    L = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
         r"$n$ & BTHM-0 meet & cap & BTHM-C meet & rescued & worst trims \\", r"\midrule"]
    for n in ("512", "1024", "2048"):
        R2 = {r: by for r, by in R.items() if r != "monotone"}
        tot = sum(c["tapes"] for by in R2.values() for k, c in by.items() if k == n)
        mt = sum(c["b0_meet"] for by in R2.values() for k, c in by.items() if k == n)
        cy = sum(c["b0_cyc"] for by in R2.values() for k, c in by.items() if k == n)
        mc = sum(c["bC_meet"] for by in R2.values() for k, c in by.items() if k == n)
        rs = sum(c["bC_rescue"] for by in R2.values() for k, c in by.items() if k == n)
        wt = max([c["bC_trim_max"] for by in R.values() for k, c in by.items()
                  if k == n and c["bC_trim_max"] is not None] or [0])
        L.append("%s & %s & %s & %s & %s & %s \\\\" % (
            n, f(mt / tot), f(cy / tot), f(mc / tot), f(rs / cy), "{:,}".format(int(wt))))
    L += [r"\bottomrule", r"\end{tabular}"]
    io.open(os.path.join(HERE, "_agg_table.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


if __name__ == "__main__":
    rep_table(); agg_table()
    print("wrote _rep_table.tex and _agg_table.tex")
