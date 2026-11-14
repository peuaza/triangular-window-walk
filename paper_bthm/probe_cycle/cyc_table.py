"""Generate the BTHM-C table from probe_cycle/cyc_cells.json.

One row per rule at n = 512: BTHM-0 meeting and cap fractions, the cycle's length and width, the BTHM-C
meeting fraction, the rescued share of the cap-diers, the median trims, and the two read medians.
The text quotes the aggregates for all three sizes; this table is the n = 512 slice because it is the one
every rule reached.
"""
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "cyc_cells.json"), encoding="utf-8"))
ORDER = ["canonical", "runs", "dense2", "peak", "trip", "oddrun", "pattern", "LAS", "MAXP", "KCHG",
         "BAL>m/2", "BAL>m/4", "NBLK=1", "NBLK>2", "NBLK>3", "AMAX", "AMIN", "longrun", "pal", "hashmod",
         "monotone"]


def fmt(x, nd=2):
    if x is None:
        return "--"
    return ("%%.%df" % nd) % x


def tex_rule(r):
    return r"\texttt{%s}" % r.replace(">", r"$>$").replace("=", r"$=$").replace("_", r"\_")


def main():
    lines = []
    lines.append(r"\begin{tabular}{@{}lrrrrrrrrr@{}}")
    lines.append(r"\toprule")
    lines.append(r"rule & \multicolumn{2}{c}{BTHM-0} & \multicolumn{2}{c}{cycle} & \multicolumn{3}{c}{BTHM-C}"
                 r" & \multicolumn{2}{c}{reads (median)} \\")
    lines.append(r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-8}\cmidrule(lr){9-10}")
    lines.append(r" & meet & cap & len & width & meet & rescued & trims & 0 & C \\")
    lines.append(r"\midrule")
    for r in ORDER:
        c = D["rules"].get(r, {}).get("512")
        if c is None:
            continue
        meet0 = fmt(c["b0_meet_frac"], 3)
        cap = fmt(c["b0_cap_frac"], 3)
        clen = fmt(c["cyc_len_med"], 0) if c["b0_cyc"] else "--"
        cw = fmt(c["cyc_width_med"], 0) if c["b0_cyc"] else "--"
        meetc = fmt(c["bC_meet_frac"], 3)
        res = fmt(c["bC_rescue_frac"], 3) if c["bC_rescue"] else "--"
        tr = fmt(c["bC_trim_med"], 0) if c["bC_rescue"] else "--"
        r0 = fmt(c["reads0_med_meet"], 0) if c["reads0_med_meet"] else "--"
        rc = fmt(c["readsC_med_all"], 0)
        lines.append(" & ".join([tex_rule(r), meet0, cap, clen, cw, meetc, res, tr, r0, rc]) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    out = os.path.join(HERE, "_cyc_table.tex")
    io.open(out, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("wrote", out, "rows", len(ORDER))


if __name__ == "__main__":
    main()
