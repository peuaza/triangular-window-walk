"""Consistency checker for the BTHM-C paper.

Four checks, in the order of the errors this project has actually made:
  A  presence   - every claim in claims_c.json appears in the paper or in an included table
  B  bans       - wording that has already been wrong once
  C  counts     - the aggregate's convention is stated
  D  refs       - every \\ref has a \\label and every \\label is used
Exit code 1 on any hit.
"""
import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.join(HERE, "bthm_c.tex")
TABLE_DIR = os.path.join(HERE, "probe_cycle")
INCLUDED = ["_cyc_table.tex", "_agg_table.tex", "_rep_table.tex"]


def main():
    text = io.open(PAPER, encoding="utf-8").read()
    for t in INCLUDED:
        p = os.path.join(TABLE_DIR, t)
        if os.path.exists(p):
            text += "\n" + io.open(p, encoding="utf-8").read()
    claims = json.load(io.open(os.path.join(HERE, "claims_c.json"), encoding="utf-8"))
    hits = []

    for c in claims["claims"]:
        for a in c.get("appear", []):
            n = len(re.findall(a["re"], text))
            if n < a.get("min", 1):
                hits.append("A presence: %s wants %r at least %d time(s), found %d" % (c["id"], a["re"], a.get("min", 1), n))

    for b in claims.get("bans", []):
        if re.search(b["re"], text, re.I):
            hits.append("B ban: %r present (%s)" % (b["re"], b.get("why", "")))

    for k in claims.get("counts", []):
        n = len(re.findall(k["re"], text))
        if n < k.get("min", 1):
            hits.append("C count: %r found %d time(s) (%s)" % (k["re"], n, k.get("why", "")))

    refs = set(re.findall(r"\\ref\{([^}]+)\}", text))
    labels = set(re.findall(r"\\label\{([^}]+)\}", text))
    for r in sorted(refs - labels):
        hits.append("D ref: \\ref{%s} has no \\label" % r)
    for l in sorted(labels - refs):
        if not l.startswith("sec:") and not l.startswith("tab:"):
            hits.append("D label: \\label{%s} never referenced" % l)

    print("BTHM-C check: %d claims, %d bans, %d counts, %d refs, %d labels"
          % (len(claims["claims"]), len(claims.get("bans", [])), len(claims.get("counts", [])), len(refs), len(labels)))
    if not hits:
        print("NO HITS")
        return 0
    for h in hits:
        print("  " + h)
    print("%d hit(s)" % len(hits))
    return 1


if __name__ == "__main__":
    sys.exit(main())
