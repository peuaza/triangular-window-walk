#!/usr/bin/env bash
# Build the merged paper, but only after its own consistency check passes.
set -euo pipefail
cd "$(dirname "$0")"

PY="C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe"
export PATH="/c/Users/shini/AppData/Local/Programs/MiKTeX/miktex/bin/x64:$PATH"

echo "== the merged source"
"$PY" _mk_ck.py

echo "== tables"
"$PY" probe_trim/cut_tables.py
"$PY" probe_colour_k/k_table.py
"$PY" probe_colour_k/k_fig.py

echo "== consistency check"
"$PY" check_ck.py

echo "== pdflatex"
for i in 1 2; do
  pdflatex -interaction=nonstopmode bthm_ck.tex >/dev/null 2>&1 || true
done
grep "Output written" bthm_ck.log
errs=$(grep -c "^!" bthm_ck.log || true)
undef=$(grep -c "undefined" bthm_ck.log || true)
echo "errors=$errs undefined=$undef"
test "$errs" -eq 0
test "$undef" -eq 0
