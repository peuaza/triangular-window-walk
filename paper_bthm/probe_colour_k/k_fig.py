"""The k-sweep figure: the meeting fraction against the number of colours the rule may read.

Nine lines, one per (rule, n). The point of the figure is where the gain sits: the big step is k = 1 -> 2 for
runs and peak, the curve flattens from k = 3, and no line is flat between 1 and 2 -- that is what "front-loaded
and decelerating, not saturated" looks like.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "k_cells.json"), encoding="utf-8"))["cells"]
RULES = ["pattern", "runs", "peak"]
NS = ["512", "1024", "2048"]
STYLE = {  # rule -> (colour, marker)
    "pattern": ("#4a6fa5", "o"),
    "runs": ("#b5651d", "s"),
    "peak": ("#5b8c5a", "^"),
}
LS = {"512": ":", "1024": "--", "2048": "-"}

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(6.6, 3.0), gridspec_kw={"width_ratios": [1.35, 1.0]})
for r in RULES:
    for n in NS:
        y = [D[r][str(k)][n]["meet_frac"] for k in (1, 2, 3, 4)]
        ax.plot([1, 2, 3, 4], y, color=STYLE[r][0], marker=STYLE[r][1], ls=LS[n], ms=3.4, lw=1.1,
                label="%s, %s" % (r, n))
ax.set_xticks([1, 2, 3, 4])
ax.set_xlabel("colours the rule may read ($k$)", fontsize=8)
ax.set_ylabel("meeting fraction", fontsize=8)
ax.tick_params(labelsize=8)
ax.set_title("the count buys meetings", fontsize=9)
ax.legend(fontsize=6.2, frameon=False, ncol=1, handlelength=1.6, labelspacing=0.25, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

# right panel: the step gain, as tapes of 300, for each (rule, n)
labels, g12, g23, g34 = [], [], [], []
for r in RULES:
    for n in NS:
        v = [D[r][str(k)][n]["meet"] for k in (1, 2, 3, 4)]
        labels.append("%s %s" % (r, n))
        g12.append(v[1] - v[0])
        g23.append(v[2] - v[1])
        g34.append(v[3] - v[2])
y = range(len(labels))
ax2.barh([i + 0.26 for i in y], g12, height=0.24, color="#4a6fa5", label="$1\\to2$")
ax2.barh([i + 0.00 for i in y], g23, height=0.24, color="#8fa8c8", label="$2\\to3$")
ax2.barh([i - 0.26 for i in y], g34, height=0.24, color="#d9d9d9", label="$3\\to4$")
ax2.set_yticks(list(y))
ax2.set_yticklabels(labels, fontsize=7)
ax2.set_xlabel("step gain, tapes of 300", fontsize=8)
ax2.tick_params(labelsize=7)
ax2.set_title("front-loaded, then flat", fontsize=9)
ax2.legend(fontsize=6.2, frameon=False, handlelength=1.2, labelspacing=0.25)
ax2.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "figs", "k_meet." + ext), dpi=200)
print("wrote figs/k_meet.pdf and .png")
