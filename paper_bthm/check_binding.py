"""Value binding: every headline number must sit in a clause whose level and direction match its source.

check_paper.py only asks whether a string occurs somewhere.  This asks a second question of the same
clause: is the number described at the level it was obtained at, and in the right direction?  A
machine-level number described as a rule-level one, or an upper bound described as a lower one, is the
error class that has recurred most often in this paper, and neither is visible to a presence check.

Scope matters.  Scanning whole paragraphs produced 23 hits, almost all false, because a paragraph may
discuss lower bounds for a different object.  Direction is therefore checked in a window of WINDOW
characters around the claim, level in the sentence containing it, and a level hit is a warning unless the
sentence also names the machine.

claims.json entries may carry "level": "rule"|"machine" and "direction": "upper"|"lower"|"exact".
Exit code 1 on any direction hit.
"""
import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WINDOW = 150

DIRECTION_FORBID = {
    "upper": [r"lower\s+bound", r"never\s+a\s+lower"],
    "lower": [r"upper\s+bound"],
    "exact": [r"upper\s+bound", r"lower\s+bound"],
}
LEVEL_FORBID = [r"at\s+the\s+level\s+of\s+the\s+rule", r"rule-level"]
MACHINE_WORDS = [r"\bmachine", r"\bdevice\b", r"\bM\^"]


def sentences(text):
    for m in re.finditer(r"[^.\n]*(?:\.[^.\n]*)*\.", text):
        yield m.start(), m.group(0)
    tail = re.search(r"[^.\n]+$", text)
    if tail:
        yield tail.start(), tail.group(0)


def main():
    text = io.open(os.path.join(HERE, "bthm_v3.tex"), encoding="utf-8").read()
    claims = json.load(io.open(os.path.join(HERE, "claims.json"), encoding="utf-8"))
    dirhits, lvlwarn = [], []
    checked = 0

    def walk(o):
        if isinstance(o, dict):
            if o.get("id") and (o.get("level") or o.get("direction")):
                yield o
            for v in o.values():
                yield from walk(v)
        elif isinstance(o, list):
            for v in o:
                yield from walk(v)

    for c in walk(claims):
        pats = [a.get("re") for a in c.get("appear", []) if a.get("re")]
        if not pats:
            continue
        checked += 1
        for p in pats:
            for m in re.finditer(p, text):
                lo, hi = max(0, m.start() - WINDOW), min(len(text), m.end() + WINDOW)
                win = text[lo:hi]
                line = text[:m.start()].count("\n") + 1
                d = c.get("direction")
                if d:
                    sent_d = None
                    for off, s in sentences(text):
                        if off <= m.start() < off + len(s) + 3:
                            sent_d = s
                            break
                    if sent_d:
                        for rx in DIRECTION_FORBID.get(d, []):
                            if re.search(rx, sent_d, re.I):
                                dirhits.append((c["id"], d, rx, line, sent_d.replace("\n", " ")[:120]))
                if c.get("level") == "machine":
                    sent = None
                    for off, s in sentences(text):
                        if off <= m.start() < off + len(s) + 3:
                            sent = s
                            break
                    if sent and any(re.search(rx, sent) for rx in LEVEL_FORBID):
                        if not any(re.search(w, sent, re.I) for w in MACHINE_WORDS):
                            lvlwarn.append((c["id"], line, sent.replace("\n", " ")[:120]))

    print("value-binding: %d claims carry a level/direction" % checked)
    for w in lvlwarn:
        print("[warn level] %s line %d: %s" % w)
    for h in dirhits:
        print("[direction] %s expects %s, clause says %r (line %d): %s" % h)
    print("%d direction hit(s), %d level warning(s)" % (len(dirhits), len(lvlwarn)))
    return 1 if dirhits else 0


if __name__ == "__main__":
    sys.exit(main())
