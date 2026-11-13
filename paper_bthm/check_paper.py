"""Data-driven consistency check for bthm_v3.tex + _appendix_evidence.tex.

Every headline number lives in claims.json with the probe that generated it; the spelling bans and the
catalogue census live there too.  This file contains no facts of its own -- only the checks:

  A  banned vocabulary (renamed classes, retired actions, refuted conjectures)
  B  census: every catalogue row present, no retired class tag
  C  each claim's pattern present the required number of times, its stale alternative absent
  D  every \\ref has a \\label

Exit 1 on any hit.  Run before pdflatex: a paper that fails this is not a paper.

Usage: python check_paper.py [-v]
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = json.load(open(os.path.join(HERE, 'claims.json'), encoding='utf-8'))
TEXT = {n: open(os.path.join(HERE, n), encoding='utf-8').read() for n in CFG['paper']}
ALL = "\n".join(TEXT.values())
VERBOSE = '-v' in sys.argv

hits = []
def hit(kind, where, what, why):
    hits.append((kind, where, what, why))

def each_line():
    for name, text in TEXT.items():
        for i, ln in enumerate(text.split('\n'), 1):
            yield name, i, ln

# ---------------------------------------------------------------- A. vocabulary
for b in CFG['banned']:
    for name, i, ln in each_line():
        if b['term'] in ln:
            if b.get('unless') and b['unless'] in ln:
                continue
            hit('A', f'{name}:{i}', b['term'], b['why'])

# the action is trim; 'cut' survives only in the companion-paper paragraph and the statistical sense
CUT_OK = ('companion paper', 'It is not the', '10\\times$-median', 'the $100\\,n$', 'median cut', 'statistical')
for name, i, ln in each_line():
    if re.search(r'\bcut\b', ln) and not any(a in ln for a in CUT_OK):
        hit('A', f'{name}:{i}', 'cut', "the action is trim; cut is the companion paper's")

# no retired class tag in the catalogue
for name, i, ln in each_line():
    for tag in ('WS', 'SS', 'QL'):
        if re.search(r'&\s*%s\s*\\\\' % tag, ln):
            hit('A', f'{name}:{i}', f'class tag {tag}', 'only FL / SL / SU exist')

# ---------------------------------------------------------------- B. census
APP = TEXT.get('_appendix_evidence.tex', '')
def _has(text, name):
    """the catalogue's row for `name`, ignoring LaTeX escapes"""
    return ('texttt{%s}' % name) in text.replace('\\', '')

for cls, rows in CFG['census'].items():
    for r in rows:
        if not _has(APP, r):
            hit('B', 'appendix', r, f'{cls} member missing from the catalogue')

# ---------------------------------------------------------------- C. claims
for c in CFG['claims']:
    for a in c.get('appear', []):
        n = len(re.findall(a['re'], ALL))
        if n < a['min']:
            hit('C', c['id'], f"{c['value']}: found {n}, need {a['min']}", c['source'])
    for f in c.get('forbid', []):
        m = re.search(f['re'], ALL)
        if m:
            if f.get('unless') and f['unless'] in ALL[max(0, m.start()-120):m.end()+120]:
                continue
            hit('C', c['id'], f"stale alternative present ({f['re']})", f['why'])

# ---------------------------------------------------------------- D. refs
labels = set(re.findall(r'\\label\{([^}]+)\}', ALL))
refs = set(re.findall(r'\\ref\{([^}]+)\}', ALL))
for r in sorted(refs - labels):
    hit('D', 'both', r, 'refers to a label that does not exist')
unused = sorted(labels - refs)

# ---------------------------------------------------------------- report
print(f"claims {len(CFG['claims'])}, banned {len(CFG['banned'])}, "
      f"labels {len(labels)}, refs {len(refs)}, unused labels {len(unused)}")
if VERBOSE and unused:
    print("  unused labels:", ', '.join(unused))
if not hits:
    print("NO HITS")
    sys.exit(0)
for k, w, what, why in hits:
    print(f"[{k}] {w}: {what}  -- {why}")
print(f"\n{len(hits)} hits")
sys.exit(1)
