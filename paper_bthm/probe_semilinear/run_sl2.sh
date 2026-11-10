#!/bin/sh
# Everything in SEMILINEAR.md, in order.  ~2 min for the table, ~12 min for the content
# probes at 2048 (they run 700 walks of 2048 steps and 5 erasure sweeps).
set -e
cd "$(dirname "$0")/.."
PY=${PY:-C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe}

$PY probe_semilinear/sl2_calib.py 512 30              # what pcr_lib charges, per rule
$PY probe_semilinear/sl2_gate.py 128,256,512,1024 40 | tee probe_semilinear/sl2_gate.log
$PY probe_semilinear/sl2_measure.py | tee probe_semilinear/sl2_measure.log
$PY probe_semilinear/sl2_table.py   | tee probe_semilinear/sl2_table.log
$PY probe_semilinear/sl2_controls.py
$PY probe_semilinear/sl2_content.py 512 40 6
$PY probe_semilinear/sl2_content.py 2048 30 4 | tee probe_semilinear/sl2_content_2048.log
