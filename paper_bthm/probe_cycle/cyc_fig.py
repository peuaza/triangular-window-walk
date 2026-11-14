"""The cycle-width figure: where in its window the walk repeats, per rule, at n = 512.

Horizontal bars of (median width at the repeat)/n, sorted, with the median period length printed beside each
bar. The point of the figure is that a cycle is not always a small-window affair: five rules repeat only
after the window has grown to a third of the tape or more, which is why the repair cascade is expensive.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "cyc_cells.json"), encoding="utf-8"))
R = D["rules"]
rows = []
for r, by in R.items():
    c = by.get("512")
    if not c or not c.get("b0_cyc"):
        continue
    rows.append((r, c["cyc_width_med"] / 512.0, c["cyc_len_med"], c["b0_cyc_frac"]))
rows.sort(key=lambda x: x[1])

fig, ax = plt.subplots(figsize=(6.2, 4.6))
y = range(len(rows))
ax.barh(list(y), [x[1] for x in rows], color="#4a6fa5", height=0.68)
for i, (r, w, L, cap) in enumerate(rows):
    ax.text(w + 0.008, i, "period %d, cap %.2f" % (L, cap), va="center", fontsize=7, color="#333333")
ax.set_yticks(list(y))
ax.set_yticklabels([x[0].replace(">", "$>$").replace("=", "$=$") for x in rows], fontsize=8)
ax.set_xlabel("median window width at the repeat, as a fraction of the tape ($n=512$)", fontsize=8)
ax.set_xlim(0, 0.62)
ax.tick_params(axis="x", labelsize=8)
ax.set_title("Where the walk repeats itself", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "figs", "cyc_width." + ext), dpi=200)
print("wrote figs/cyc_width.pdf and .png,", len(rows), "rules")
