#!/usr/bin/env bash
# Build the paper, but only after the consistency check passes.
# A paper that fails check_paper.py is not built: consistency is a precondition, not a review step.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-python3}"
export PATH="/c/Users/shini/AppData/Local/Programs/MiKTeX/miktex/bin/x64:$PATH"

echo "== consistency check"
"$PY" check_paper.py

echo "== pdflatex"
for i in 1 2; do
  pdflatex -interaction=nonstopmode bthm_jalc_v2.tex >/dev/null 2>&1 || true
done
grep "Output written" bthm_jalc_v2.log
errs=$(grep -c "^!" bthm_jalc_v2.log || true)
undef=$(grep -c "undefined" bthm_jalc_v2.log || true)
echo "errors=$errs undefined=$undef"
test "$errs" -eq 0
test "$undef" -eq 0
