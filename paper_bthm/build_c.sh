#!/usr/bin/env bash
# Build the BTHM-C paper, but only after its own consistency check passes.
set -euo pipefail
cd "$(dirname "$0")"

PY="C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe"
export PATH="/c/Users/shini/AppData/Local/Programs/MiKTeX/miktex/bin/x64:$PATH"

echo "== tables"
"$PY" probe_cycle/cyc_table.py
"$PY" probe_cycle/cyc_tables2.py

echo "== consistency check"
"$PY" check_c.py

echo "== pdflatex"
for i in 1 2; do
  pdflatex -interaction=nonstopmode bthm_c.tex >/dev/null 2>&1 || true
done
grep "Output written" bthm_c.log
errs=$(grep -c "^!" bthm_c.log || true)
undef=$(grep -c "undefined" bthm_c.log || true)
echo "errors=$errs undefined=$undef"
test "$errs" -eq 0
test "$undef" -eq 0
