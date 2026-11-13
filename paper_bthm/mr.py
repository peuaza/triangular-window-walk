"""Meeting rate: of the tapes drawn from the ensemble, how many does the walk meet?

The paper reports costs on meeting tapes; this records the draw count and the met count so the
conditioning is visible.  Reads nothing but tr_lib.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, "probe_triple_nat"))
import tr_lib as T  # noqa: E402

RULES = ["monotone", "pattern", "runs", "peak", "pal", "trip", "NBLK=1"]
out = {}
for n in (512, 1024):
    for rule in RULES:
        met = 0
        tot = 0
        for i, a, s, X in T.tapes(n):
            r = T.run(rule, a, s)
            tot += 1
            if r.get("delta") is not None:
                met += 1
        out["%s@%d" % (rule, n)] = {"met": met, "drawn": tot, "rate": met / tot}
        print(rule, n, met, tot, round(met / tot, 4), flush=True)
json.dump(out, open(os.path.join(HERE, "meet_rate.json"), "w"), indent=1)
print("wrote meet_rate.json")
