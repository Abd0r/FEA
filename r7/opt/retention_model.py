"""Room-temperature retention of a single-electron DB bit (1 = DB-, 0 = DB0) and the effect of each optimisation.
Every number printed comes from here. Physics (per DB, T = 300 K):
  1->0  electron emission to the CB:  e_n = nu * exp(-(dE + Ecap_n)/kT)          [detailed balance, multiphonon]
  0->1  (a) capture of a CB electron: c_n = sig_inf*exp(-Ecap_n/kT) * v_th * n_s
        (b) hole emission to the VB:  e_p = nu * exp(-(Eg - dE + Ecap_p)/kT)
  both  (c) tunnel exchange with the nearest electron reservoir at distance d: G_t = nu_t * exp(-2 kappa d)
  Ecap (opt 1) = Marcus/harmonic configuration-coordinate capture barrier (lam - D)^2/(4 lam), D = level depth.
Optimisations:
  1 lattice relaxation (lam > 0)   2 depleted/intrinsic layer (n_s small)   3 tunnel spacer (d large)
  4 level tuning by a fixed neighbour charge: buried P+ (deeper) or like charge (shallower, DFT tilt 0.081 eV at 3.8 A).
    Modelled as freedom to place the level within +-0.30 eV of mid-gap (14.4/(11.7*0.4 nm) ~ 0.3 eV); without it the level is mid-gap.
  5 dual-rail pair: a lost electron becomes a detectable erasure -> code corrects 2t erasures instead of t errors
  6 whole-Zone FZC refresh + BCH: interval and power for 8 GB."""
import math, json, itertools
from math import comb
kT = 8.617e-5*300; Eg = 1.12; vth = 2.3e7; Nc = 2.8e19          # cm/s, cm^-3 (Si, 300 K)
HBAR, ME, Q = 1.0546e-34, 9.109e-31, 1.602e-19

def kappa(U, m=0.26):                                             # 1/nm, evanescent decay under barrier U (eV)
    return math.sqrt(2*m*ME*U*Q)/HBAR*1e-9

def ecap(lam, D):
    return 0.0 if lam <= 0 else (lam - D)**2/(4*lam)

def rates(dE=0.56, lam=0.0, n_s=1e18, d_nm=1.0, sig=1e-15, nu_fixed=1e13, dq=0.0):
    dE = dE + dq
    nu = nu_fixed if lam <= 0 else sig*vth*Nc                      # baseline keeps the 1e13 phonon bound
    en = nu*math.exp(-(dE + ecap(lam, dE))/kT)
    cn = sig*math.exp(-ecap(lam, dE)/kT)*vth*n_s
    ep = nu*math.exp(-(Eg - dE + ecap(lam, Eg - dE))/kT)
    gt = 1e13*math.exp(-2*kappa(min(dE, Eg - dE))*d_nm)
    return dict(r10=en + gt, r01=cn + ep + gt)

def tau(**kw):
    r = rates(**kw); return 1/max(r['r10'], r['r01'])

def refresh(tau_s, dual=False, bits=8e9*8, t=4, n=92, E_aJ=(1, 35)):
    """Largest interval T with < 1 uncorrectable word per hour in 8 GB; dual-rail doubles correctable count."""
    tt = 2*t if dual else t; words = bits/64
    lo, hi = 1e-15, tau_s
    for _ in range(200):
        T = math.sqrt(lo*hi); p = min(1.0, T/tau_s)
        fail = comb(n, tt + 1)*p**(tt + 1)*words/T*3600
        lo, hi = (T, hi) if fail < 1 else (lo, T)
    area = 2 if dual else 1
    return T, [e*1e-18*bits*(n/64)*area/T for e in E_aJ]

if __name__ == '__main__':
    import sys
    lam = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5        # replaced by the DFT value when cc.json lands
    base = dict(dE=0.56, lam=0.0, n_s=1e18, d_nm=1.0, dq=0.0)
    opts = {1: dict(lam=lam), 2: dict(n_s=1e10), 3: dict(d_nm=10.0), 4: dict()}
    names = {1: 'lattice relaxation', 2: 'depleted layer', 3: 'tunnel spacer 10 nm', 4: 'level tuning by neighbour charge',
             5: 'dual-rail erasures', 6: 'FZC refresh'}
    rows = []
    for k in range(0, 5):
        for combo in itertools.combinations([1, 2, 3, 4], k):
            p = dict(base)
            for c in combo: p.update(opts[c])
            levels = [x/100 for x in range(26, 87)] if 4 in combo else [0.56]
            best = max((tau(**{**p, 'dE': x}), x) for x in levels)
            for dual in (False, True):
                T, P = refresh(best[0], dual)
                rows.append(dict(opts=list(combo) + ([5] if dual else []), dE=best[1], tau_s=best[0], refresh_s=T, W_1aJ=P[0], W_35aJ=P[1]))
    print(f'lambda = {lam} eV  (opt 1 input)')
    print(f"{'optimisations':28s} {'level':>5s} {'bit lifetime':>13s} {'refresh every':>14s} {'W @1aJ':>9s} {'W @35aJ':>9s}")
    for r in rows:
        lab = '+'.join(map(str, r['opts'])) or 'none (baseline)'
        print(f"{lab:28s} {r['dE']:5.2f} {r['tau_s']:13.3g} {r['refresh_s']:14.3g} {r['W_1aJ']:9.3g} {r['W_35aJ']:9.3g}")
    json.dump(dict(lam=lam, rows=rows, names=names), open(f'retention_lam{lam}.json', 'w'), indent=1)
