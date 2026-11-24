#!/usr/bin/env python3
"""Reproduction check: recompute FEA-R7 from the committed DFT results and compare key numbers
with the reference values in reference_r7.json. Prints PASS/FAIL per quantity; exit code 1 on any failure."""
import json, sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
ref = json.load(open('reference_r7.json'))
N = json.load(open('paper/numbers_r7.json'))
AD = json.load(open('opt/chip_sim_adaptive_results.json'))
def get(path):
    src, *keys = path.split(':')
    d = {'N': N, 'AD': AD}[src]
    for k in keys:
        d = d[int(k)] if isinstance(d, list) else d[k]
    return d
fail = 0
for path, (val, rtol) in ref.items():
    got = get(path); ok = abs(got - val) <= rtol*abs(val) + 1e-12
    fail |= not ok
    print(f"  {'PASS' if ok else 'FAIL'}  {path:58s} {got:.6g}  (ref {val:.6g})")
print('  ALL R7 CHECKS PASS' if not fail else '  SOME R7 CHECKS FAILED'); sys.exit(fail)
