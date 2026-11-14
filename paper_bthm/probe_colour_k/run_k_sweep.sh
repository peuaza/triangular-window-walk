#!/bin/bash
# The k-sweep grid: 3 rules x 3 sizes x k in {1,2,3,4}, 300 tapes each,
# 8 workers (16 cores; numba is single-threaded per cell).
cd "$(dirname "$0")"
PY=../probe_cycle/.venv/Scripts/python.exe
LOG=k_run.log
: > "$LOG"
run() {
  echo "START $2 n=$1 k=$3 $(date +%H:%M:%S)" >> "$LOG"
  "$PY" k_run.py "$1" "$2" "$3" 300 >> "$LOG" 2>&1
  echo "END   $2 n=$1 k=$3 code=$? $(date +%H:%M:%S)" >> "$LOG"
}
for n in 512 1024 2048; do
  for rule in pattern runs peak; do
    for k in 4 2 3 1; do
      run "$n" "$rule" "$k" &
      while [ "$(jobs -r | wc -l)" -ge 8 ]; do wait -n; done
    done
  done
done
wait
touch DONE
echo "ALL DONE $(date +%H:%M:%S)" >> "$LOG"
