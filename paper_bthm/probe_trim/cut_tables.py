"""Tables for BTHM-C from the cut measurement (probe_trim/cut_cells.json).

_trim_table.tex    per rule at n = 512: BTHM-0 meet/cap, cycle length/width, cut-repair meet/empty/rescued,
                  trims median/max, reads medians
_trimagg_table.tex the aggregate at the three sizes over the twenty rules that were run
_trimrep_table.tex repair semantics: window retraction (probe_cycle) against cut (probe_trim), n = 512
"""
import io
import json
import os

CUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "probe_trim", "cut_cells.json")
OLD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "probe_cycle", "cyc_cells.json")
D = json.load(open(CUT, encoding="utf-8"))
O = json.load(open(OLD, encoding="utf-8"))
R = D["rules"]
ORDER = ["canonical", "runs", "dense2", "peak", "trip", "oddrun", "pattern", "LAS", "MAXP", "KCHG",
         "BAL>m/2", "BAL>m/4", "NBLK=1", "NBLK>2", "NBLK>3", "AMAX", "AMIN", "longrun", "pal", "hashmod",
         "monotone"]


def f(x, nd=3):
    return "--" if x is None else ("%%.%df" % nd) % x


def i0(x):
    return "--" if x is None else "%d" % round(x)


def nm(r):
    return r"\texttt{%s}" % r.replace(">", r"$>$").replace("=", r"$=$").replace("_", r"\_")


def per_rule():
    L = [r"\begin{tabular}{@{}lrrrrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{BTHM-0} & \multicolumn{2}{c}{cycle} & \multicolumn{3}{c}{BTHM-C (cut)}"
         r" & \multicolumn{2}{c}{reads} \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-8}\cmidrule(lr){9-10}",
         r"rule & meet & cap & len & width & meet & rescued & trims & 0 & C \\", r"\midrule"]
    for r in ORDER:
        c = R.get(r, {}).get("512")
        if not c:
            continue
        L.append(" & ".join([nm(r), f(c["b0_meet_frac"]), f(c["b0_cyc_frac"]),
                             i0(c["cyc_len_med"]) if c["b0_cyc"] else "--",
                             i0(c["cyc_width_med"]) if c["b0_cyc"] else "--",
                             f(c["bC_meet_frac"]), f(c["bC_rescue_frac"]) if c["bC_rescue"] else "--",
                             i0(c["bC_trim_med"]) if c["bC_rescue"] else "--",
                             i0(c["reads0_med_meet"]) if c["reads0_med_meet"] else "--",
                             i0(c["readsC_med_all"])]) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    io.open("_trim_table.tex", "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


def agg():
    L = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
         r"$n$ & BTHM-0 meet & cap & BTHM-C meet & rescued & max trims \\", r"\midrule"]
    for n in ("512", "1024", "2048"):
        keys = ("tapes", "b0_meet", "b0_cyc", "bC_meet", "bC_rescue")
        s = {k: 0 for k in keys}
        for r, by in R.items():
            if r == "monotone":
                continue
            c = by.get(n)
            if not c:
                continue
            for k in keys:
                s[k] += c.get(k) or 0
        tm = max([c["bC_trim_max"] for r, by in R.items() if r != "monotone"
                  for k, c in by.items() if k == n and c.get("bC_trim_max") is not None] or [0])
        L.append("%s & %s & %s & %s & %s & %s \\\\" % (
            n, f(s["b0_meet"] / s["tapes"]), f(s["b0_cyc"] / s["tapes"]), f(s["bC_meet"] / s["tapes"]),
            f(s["bC_rescue"] / s["b0_cyc"]), "{:,}".format(int(tm))))
    L += [r"\bottomrule", r"\end{tabular}"]
    io.open("_trimagg_table.tex", "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


def repair():
    L = [r"\begin{tabular}{@{}lrrr@{}}", r"\toprule",
         r"rule & retraction & cut & $\Delta$ \\", r"\midrule"]
    for r in ORDER:
        if r == "monotone":
            continue
        a = O["rules"].get(r, {}).get("512", {}).get("bC_meet_frac")
        b = R.get(r, {}).get("512", {}).get("bC_meet_frac")
        if a is None or b is None:
            continue
        L.append("%s & %s & %s & %+.3f \\\\" % (nm(r), f(a), f(b), b - a))
    L += [r"\bottomrule", r"\end{tabular}"]
    io.open("_trimrep_table.tex", "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)
    per_rule(); agg(); repair()
    print("wrote _trim_table.tex, _trimagg_table.tex, _trimrep_table.tex")
