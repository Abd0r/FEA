"""Whole-chip FEA simulation at room temperature with all five optimisations applied.

Chip (from R7, ../paper/numbers_r7.json 'chip' block, written by fea_r7.py). Earlier R6 version: one injection
= 1 eV = 0.16 aJ; control-line switching 35 aJ per wire operation; line RC 40 ns.
Usable working memory = raw x 64/n for the BCH code in the chip block (R7 design: n=120, t=8, d=17,
so up to 16 erasures or 8 errors per word).

Optimisations (all on):
  1 lattice-relaxation capture barrier, Marcus form, lambda from DFT (scanned until cc.json lands)
  2 depleted / intrinsic layer: free electrons at the Block = n_i(T)
  3 tunnel spacer: 10 nm to the nearest reservoir
  4 level tuning by a fixed neighbour charge: level placed at the best depth (found per lambda at 300 K, then frozen)
  5 dual-rail bits: a lost electron is a detectable erasure
  + whole-Zone FZC refresh with an adaptive interval: the FZC reads chip temperature and picks the longest
    interval that keeps the chip under ONE uncorrectable word per year.

Real-world use: a 24-hour day of mixed workloads, in 3 machines (phone, laptop, desktop), with a lumped
thermal model (T = T_amb + R_th * P, time constant tau_th). Refresh power rises with temperature, which
heats the chip: the loop is solved each step and a thermal runaway is flagged if no stable point exists.
A Monte Carlo of dual-rail words checks the analytic decoder-failure formula.

Usage: python3 fea_chip_sim.py [lambda_eV ...]
"""
import math, json, sys
import numpy as np
from math import comb

KB = 8.617e-5
_here = __import__('os').path.dirname(__import__('os').path.abspath(__file__))
CHIP = json.load(open(_here + '/../paper/numbers_r7.json'))['chip']
RAW_BITS = CHIP['raw_bits']                      # dual-rail pairs = logical raw bits
CODE_N, CODE_D = CHIP['code_n'], CHIP['code_d']
WORDS = RAW_BITS/CODE_N
USABLE_GB = WORDS*64/8e9
INJ_J = 1.602176634e-19
# Read of an electron bit = one R5 logic-mode pass. Two error types:
#   silent misread (empty Block passes):            READ_ERR
#   read disturbance (empty Block captures, '11'):   DISTURB -> a detectable erasure in dual rail
READ_ERR = CHIP['read_silent_err']
DISTURB = CHIP['read_disturb']
E_REFRESH_NOM = CHIP['E_refresh_J']
E_REFRESH_PESS = CHIP['E_refresh_R6lines_J']     # R6 control lines (3.3 V swing), no line improvement
E_GATE = CHIP['E_gate_J']
TARGET_FAILS_PER_YEAR = 1.0
YEAR = 3.156e7

# ---------- device physics (same model as retention_model.py, with temperature) ----------
def Eg(T):    return 1.170 - 4.73e-4*T*T/(T + 636)                 # Varshni, Si
def ni(T):    return 1.0e10*(T/300)**1.5*math.exp(-(Eg(T)/(2*KB*T) - Eg(300)/(2*KB*300)))
def vth(T):   return 2.3e7*math.sqrt(T/300)
def Nc(T):    return 2.8e19*(T/300)**1.5
def kappa(U): return math.sqrt(2*0.26*9.109e-31*U*1.602e-19)/1.0546e-34*1e-9
def ecap(lam, D): return (lam - D)**2/(4*lam)

def tau_bit(T, lam, dE, sig=1e-15, d_nm=10.0):
    kT = KB*T; g = Eg(T); nu = sig*vth(T)*Nc(T)
    en = nu*math.exp(-(dE + ecap(lam, dE))/kT)                          # 1 -> 0
    cn = sig*math.exp(-ecap(lam, dE)/kT)*vth(T)*ni(T)                   # 0 -> 1 from CB
    ep = nu*math.exp(-(g - dE + ecap(lam, g - dE))/kT)                  # 0 -> 1 from VB
    gt = 1e13*math.exp(-2*kappa(min(dE, g - dE))*d_nm)                  # tunnel exchange
    return 1/max(en + gt, cn + ep + gt)

def best_level(lam):
    return max((tau_bit(300, lam, x/100), x/100) for x in range(26, 87))[1]

# ---------- code: errors e and erasures s are correctable while 2e + s <= d - 1 ----------
CODE_K = CHIP.get('code_k', 64)
ML = CHIP.get('decoder') == 'ML'
def ml_fail(s, n=None, k=None):
    """P(an s-erasure pattern defeats Gaussian-elimination decoding): 0 below d, else 2^-(n-k-s)
    (random-code value; Monte Carlo on the real code agrees or lies below it, opt/bch_erasure.py)."""
    n = CODE_N if n is None else n; k = CODE_K if k is None else k
    return 0.0 if s < CODE_D else min(1.0, 2.0**(-(n - k - s)))

def word_fail(p_er, p_err, n=None, d=None, ml=None):
    """P(word lost). Bounded-distance: 2e + s >= d. With ML erasure decoding (default from the chip block when
    n, d are the chip's): error-free words with s >= d fail only with probability ml_fail(s)."""
    ml = (ML and n is None and d is None) if ml is None else ml
    n = CODE_N if n is None else n; d = CODE_D if d is None else d; emax = (d - 1)//2
    tot = 0.0
    for e in range(0, emax + 1):
        pe = comb(n, e)*p_err**e*(1 - p_err)**(n - e) if p_err > 0 else (1.0 if e == 0 else 0.0)
        if pe < 1e-300: continue
        smin = max(0, d - 2*e)
        # tail of erasures among the remaining n-e bits
        m = n - e; ps = 0.0
        for s in range(smin, m + 1):
            ps += comb(m, s)*p_er**s*(1 - p_er)**(m - s)
        tot += pe*ps
    for e in range(emax + 1, emax + 4):        # more than emax errors alone already fails
        tot += comb(n, e)*p_err**e
    if ml:                                     # error-free words: replace the s >= d tail by its ML failure share
        p0 = (1 - p_err)**n; dt = 0.0
        for s in range(d, n + 1):
            bs = comb(n, s)*p_er**s*(1 - p_er)**(n - s)
            dt += bs*(1 - ml_fail(s, n))
        tot = max(tot - p0*dt, 0.0)
    return tot

from functools import lru_cache

def pick_interval(T, lam, dE, reads_per_word_per_s):
    return _pick(round(T, 1), lam, dE, reads_per_word_per_s)

@lru_cache(maxsize=None)
def _pick(T, lam, dE, reads_per_word_per_s):
    """Failures vs interval are U-shaped: long intervals lose electrons (erasures), short ones multiply
    read errors (every refresh reads every bit). Scan a log grid; return the LONGEST interval that meets
    the target, or the best achievable interval if none does (flag via fails_per_s)."""
    tau = tau_bit(T, lam, dE)
    def fails(Tr):
        p_er = -math.expm1(-Tr/tau)
        return (word_fail(p_er + DISTURB, READ_ERR)*WORDS/Tr
                + word_fail(p_er/2 + DISTURB, READ_ERR)*WORDS*reads_per_word_per_s)
    grid = [10**(x/20) for x in range(-200, 61)]          # 1e-10 .. 1e3 s
    grid = [g for g in grid if g < 10*tau] or [tau]
    ok = [g for g in grid if fails(g)*YEAR < TARGET_FAILS_PER_YEAR]
    Tr = max(ok) if ok else min(grid, key=fails)
    return Tr, tau, fails(Tr)

# ---------- workloads (per-second rates for the whole chip) ----------
# traffic = words read+written by programs; ops = bitwise gate ops computed in-array; io = bits off-chip
WORK = {
    'sleep':        dict(traffic_GBs=0.05, ops=1e9,  io_GBs=0.0),
    'office/web':   dict(traffic_GBs=3,    ops=5e12, io_GBs=0.2),
    'video call':   dict(traffic_GBs=8,    ops=3e13, io_GBs=0.5),
    'gaming':       dict(traffic_GBs=60,   ops=5e14, io_GBs=4),
    'compile':      dict(traffic_GBs=25,   ops=2e14, io_GBs=0.5),
    'AI inference': dict(traffic_GBs=400,  ops=4e15, io_GBs=0.1),   # 7B-parameter model, bandwidth-bound
}
DAY = [('sleep', 8), ('office/web', 6), ('video call', 2), ('compile', 2), ('gaming', 3), ('AI inference', 3)]
MACHINES = {   # ambient (C), thermal resistance (K/W), thermal time constant (s), boundary periphery power (W)
    'phone':   dict(Tamb=30, Rth=8.0, tau_th=60,  P_edge=0.3, power_cap=5),
    'laptop':  dict(Tamb=35, Rth=1.2, tau_th=40,  P_edge=1.5, power_cap=45),
    'desktop': dict(Tamb=25, Rth=0.35, tau_th=30, P_edge=3.0, power_cap=250),
}
IO_J_PER_BIT = 5e-12          # boundary transistors driving off-chip links (the only CMOS)

def power(T, lam, dE, w, e_ref, m):
    reads = w['traffic_GBs']*8e9/64/WORDS                               # accesses per word per s
    Tr, tau, _ = pick_interval(T, lam, dE, reads)
    P_ref = RAW_BITS*e_ref/Tr
    P_ops = w['ops']*E_GATE
    P_acc = w['traffic_GBs']*8e9*(E_REFRESH_NOM)                        # an access costs about one refresh
    P_io = w['io_GBs']*8e9*IO_J_PER_BIT + m['P_edge']
    return P_ref + P_ops + P_acc + P_io, dict(P_refresh=P_ref, P_compute=P_ops, P_access=P_acc, P_io=P_io,
                                              refresh_s=Tr, tau_s=tau)

def steady(lam, dE, w, e_ref, m):
    """Fixed point of T = Tamb + Rth * P(T); runaway if P grows faster than 1/Rth."""
    T = m['Tamb'] + 273.15
    for _ in range(200):
        P, parts = power(T, lam, dE, w, e_ref, m)
        Tn = m['Tamb'] + 273.15 + m['Rth']*P
        if Tn > 273.15 + 150: return None
        if abs(Tn - T) < 1e-3: return Tn, P, parts
        T = 0.5*T + 0.5*Tn
    return T, P, parts

def day(lam, dE, e_ref, m, dt=10.0):
    T = m['Tamb'] + 273.15; log = []; fails = 0.0; energy = 0.0; runaway = None
    for name, hours in DAY:
        w = WORK[name]
        for _ in range(int(hours*3600/dt)):
            P, parts = power(T, lam, dE, w, e_ref, m)
            if P > 50*m['power_cap'] or T > 273.15 + 150:
                runaway = name; break
            Tss = m['Tamb'] + 273.15 + m['Rth']*P
            T += (Tss - T)*(1 - math.exp(-dt/m['tau_th']))
            energy += P*dt
            reads = w['traffic_GBs']*8e9/64/WORDS
            fails += _pick(round(T, 1), lam, dE, reads)[2]*dt
        if runaway: break
        log.append(dict(work=name, T_C=T - 273.15, P_W=P, **{k: v for k, v in parts.items()}))
    return dict(log=log, energy_Wh=energy/3600, expected_fails=fails, runaway=runaway)

def monte_carlo(p_er, words=200_000, rng=np.random.default_rng(1)):
    """Dual-rail words: each rail pair loses its electron with p_er (erasure); reads err with READ_ERR."""
    s = rng.binomial(CODE_N, min(1.0, p_er + DISTURB), size=words); e = rng.binomial(CODE_N - s, READ_ERR)
    fail = 2*e + s >= CODE_D
    if ML:
        pf = np.array([ml_fail(int(x)) for x in s])
        fail = np.where(e == 0, rng.random(words) < pf, fail)
    return float(np.mean(fail))

if __name__ == '__main__':
    lams = [float(x) for x in sys.argv[1:]] or [0.5, 1.0, 1.5]
    out = dict(chip=dict(raw_GB=RAW_BITS/8e9, usable_GB=USABLE_GB, words=WORDS, E_refresh_nom_aJ=E_REFRESH_NOM*1e18,
                         E_gate_aJ=E_GATE*1e18, read_err=READ_ERR), runs={})
    print(f"chip: {RAW_BITS/8e9:.1f} GB raw dual-rail working memory -> {USABLE_GB:.1f} GB usable after BCH; "
          f"refresh {E_REFRESH_NOM*1e18:.2f} aJ/bit (pessimistic 35), gate op {E_GATE*1e18:.2f} aJ")
    # decoder check
    for p in (0.30, 0.35, 0.40):          # high enough that the Monte Carlo sees ML failures
        print(f"  decoder check p_erasure={p}: analytic {word_fail(p + DISTURB, READ_ERR):.3e}  Monte Carlo {monte_carlo(p):.3e}")
    for lam in lams:
        dE = best_level(lam)
        print(f"\n=== lambda = {lam} eV, level tuned to {dE:.2f} eV below CB; bit lifetime "
              f"{tau_bit(300, lam, dE):.3g} s at 27 C, {tau_bit(323, lam, dE):.3g} s at 50 C, {tau_bit(358, lam, dE):.3g} s at 85 C")
        out['runs'][lam] = {}
        for e_lab, e_ref in (('nominal', E_REFRESH_NOM), ('pessimistic', E_REFRESH_PESS)):
            for mname, m in MACHINES.items():
                r = day(lam, dE, e_ref, m)
                out['runs'][lam][f'{mname}/{e_lab}'] = r
                if r['runaway']:
                    print(f"  {mname:8s} {e_lab:11s} THERMAL RUNAWAY during {r['runaway']}"); continue
                peak = max(r['log'], key=lambda x: x['P_W'])
                print(f"  {mname:8s} {e_lab:11s} day energy {r['energy_Wh']:8.2f} Wh | peak {peak['P_W']:7.2f} W "
                      f"at {peak['T_C']:5.1f} C ({peak['work']}) | refresh at peak every {peak['refresh_s']:.2g} s | "
                      f"expected uncorrectable words/day {r['expected_fails']:.2g}")
                for x in r['log']:
                    print(f"      {x['work']:13s} T {x['T_C']:5.1f} C  P {x['P_W']:8.3f} W  "
                          f"(refresh {x['P_refresh']:.3g}, compute {x['P_compute']:.3g}, access {x['P_access']:.3g}, io {x['P_io']:.3g})"
                          f"  refresh every {x['refresh_s']:.2g} s")
    json.dump(out, open('chip_sim_results.json', 'w'), indent=1, default=float)
