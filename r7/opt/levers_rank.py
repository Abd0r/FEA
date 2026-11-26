"""Rank temperature levers against the B+C design (16-erasure code, bounded-distance decoding, one-step
L1 lines) at the DFT lambda. E1 and E2 are now part of the R7 design; this script keeps the ablation.
Metric: highest ambient (C) at which an adaptive chip meets <=1 uncorrectable word/yr while gaming.
Usage: python3 levers_rank.py  -> levers_rank.json"""
import sys, math, json, os
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fea_chip_sim as B, fea_chip_sim_adaptive as A
LAM = round(json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'paper', 'numbers_r7.json')))['lambda_eV']['minus'], 3)
wf0, tb0 = B.word_fail, B.tau_bit
RTH0 = {k: v['Rth'] for k, v in B.MACHINES.items()}
S1, S3 = A.SENSE[1], A.SENSE[3]
E_NOM = B.CHIP['E_refresh_conv_J']           # B+C: injection + 2 one-step line switches at 3.5 aJ
LINE_NOM = (E_NOM - B.INJ_J)/2

def tb_wide(T, lam, dE, sig, shift):
    kT = B.KB*T; g = B.Eg(T) + shift; D = dE + shift/2; nu = sig*B.vth(T)*B.Nc(T)
    en = nu*math.exp(-(D + B.ecap(lam, D))/kT); ep = nu*math.exp(-(g - D + B.ecap(lam, g - D))/kT)
    cn = sig*math.exp(-B.ecap(lam, D)/kT)*B.vth(T)*B.ni(T)*math.exp(-shift/(2*kT))
    gt = 1e13*math.exp(-2*B.kappa(min(D, g - D))*10.0)
    return 1/max(en + gt, cn + ep + gt)

def setup(sig=1e-15, ml=False, line_div=1.0, rth=1.0, gap=0.0):
    B.word_fail = wf0; B.ML = ml
    B.tau_bit = (lambda T, lam, dE: tb0(T, lam, dE, sig=sig)) if gap == 0 else (lambda T, lam, dE: tb_wide(T, lam, dE, sig, gap))
    e = B.INJ_J + 2*LINE_NOM/line_div
    A.SENSE[1] = (S1[0], e); A.SENSE[3] = (S3[0], e + 2*B.INJ_J)
    for k in B.MACHINES: B.MACHINES[k]['Rth'] = RTH0[k]*rth
    A.pick.cache_clear()

V = [('B+C (R7 before E1, E2)', {}),
     ('E1 ML erasure decoding of the same code', dict(ml=True)),
     ('E2 4-step adiabatic line drive (lines /4)', dict(line_div=4)),
     ('E3 package thermal resistance -25%', dict(rth=0.75)),
     ('P1 capture cross-section 1e-17 cm2', dict(sig=1e-17)),
     ('P2 gap +0.4 eV (thin-SOI confinement)', dict(gap=0.4)),
     ('E1+E2', dict(ml=True, line_div=4)),
     ('E1+E2+E3', dict(ml=True, line_div=4, rth=0.75)),
     ('E1+E2+E3+P1', dict(ml=True, line_div=4, rth=0.75, sig=1e-17)),
     ('all (E1+E2+E3+P1+P2)', dict(ml=True, line_div=4, rth=0.75, sig=1e-17, gap=0.4))]
out = {'lambda_eV': LAM, 'rows': []}
for name, kw in V:
    setup(**kw)
    amb = {m: A.max_ambient(m, LAM, True) for m in ('phone', 'laptop', 'desktop')}
    out['rows'].append(dict(name=name, **kw, **amb))
    print(f"{name:44s} phone {amb['phone']:6.1f}  laptop {amb['laptop']:6.1f}  desktop {amb['desktop']:6.1f}", flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'levers_rank.json'), 'w'), indent=1)
