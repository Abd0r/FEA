#!/usr/bin/env python3
"""FEA R5 analysis. Every number quoted in FEA-R5.tex is written to numbers_r5.json.
Usage: python3 fea_r5.py [--figs]     (numpy; figures also need matplotlib)
"""
import json, math, sys
import numpy as np

kB = 8.617333262e-5; hbar_eVs = 6.582119569e-16; hbar = 1.054571817e-34
e = 1.602176634e-19; COUL = 1.4399645   # e^2/(4 pi eps0), eV nm
YEAR = 3.15576e7; T0 = 300.0

# ------------------------------------------------------------------ inputs
a, b, c = 0.384, 0.768, 0.225        # nm, Si(100)-2x1:H lattice
t = 0.020        # eV, wire hopping (input)
# cluster-to-wire coupling is DERIVED (Section 2 of the paper): every non-addressed
# cell sits at the hold detuning D_off and must scatter little enough that an
# electron crosses all CELLS_PER_WIRE cells with total loss <= LOSS_WIRE.
CELLS_PER_WIRE = 128      # one 64-bit dual-rail Word per wire
D_OFF = 1.0               # eV, hold/off detuning (= storage barrier scale)
LOSS_WIRE = 1e-3
Gamma_design = 2*D_OFF*math.sqrt(LOSS_WIRE/(CELLS_PER_WIRE - 1))   # R_off ~ (Gamma/2D)^2
tw = math.sqrt(Gamma_design*t)
tc = 0.015       # eV, intra-cluster hopping (input)
alpha_g = 0.30   # eV/V, gate lever arm (input)
nu = 15.6e12     # Hz, Si zone-centre optical phonon
nu_lo, nu_hi = 1e12, 15.6e12
hw_ph = 0.0645   # eV, Si optical phonon energy
area = 0.5       # cm^2 array
scrub = 64e-3    # s
N = {}

# ------------------------------------------------------------------ 1. side-coupled lossy resonance on a 1D wire
def scat(E, E0, gam, tw=tw, t=t):
    """Wire site 0 carries the self-energy of a side-coupled level E0 - i*gam.
    Returns T, R, A for an electron at energy E inside the band |E|<2t."""
    E = np.asarray(E, dtype=float)
    v = np.sqrt(np.clip(4*t**2 - E**2, 1e-30, None))      # 2 t sin k
    V = tw**2/(E - E0 + 1j*gam)
    tamp = 1j*v/(1j*v - V)
    r = tamp - 1
    Tt, Rr = np.abs(tamp)**2, np.abs(r)**2
    return Tt, Rr, 1 - Tt - Rr

Gamma_el = tw**2/t              # elastic width at band centre (= 2 tw^2 / v, v = 2t)
kappa = Gamma_el/2
N['Gamma_el_meV'] = Gamma_el*1e3
N['kappa_meV'] = kappa*1e3
N['tau_elastic_fs'] = hbar_eVs/Gamma_el*1e15
gams = np.logspace(-12, 0, 4000)*kappa
Tg, Rg, Ag = scat(0.0, 0.0, gams)
i = int(np.argmax(Ag))
N['A_max'] = float(Ag[i]); N['gamma_at_Amax_over_kappa'] = float(gams[i]/kappa)
T1, R1, A1 = scat(0.0, 0.0, kappa)
N['critical'] = dict(T=float(T1), R=float(R1), A=float(A1))
# realistic phonon-emission relaxation times
for tau_ps in (0.1, 1.0, 10.0):
    g = hbar_eVs/(tau_ps*1e-12)/2      # gamma = hbar/(2 tau)
    Tt, Rr, Aa = scat(0.0, 0.0, g)
    N[f'write_tau{tau_ps}ps'] = dict(gamma_meV=g*1e3, T=float(Tt), R=float(Rr), A=float(Aa))
# logic mode: requirement on gamma for T<=1e-9 (gate error) and A<=1e-6 (read disturb)
def need_gamma(target, which):
    for g in gams:
        Tt, Rr, Aa = scat(0.0, 0.0, g)
        val = Tt if which == 'T' else Aa
        if val > target:
            return float(prev)
        prev = g
    return float(gams[-1])
gT = need_gamma(1e-9, 'T'); gA = need_gamma(1e-6, 'A')
N['logic_gammaT_over_kappa'] = gT/kappa
N['logic_gammaA_over_kappa'] = gA/kappa
N['logic_relax_time_s'] = hbar_eVs/(2*gA)
Tl, Rl, Al = scat(0.0, 0.0, gA)
N['logic_mode'] = dict(T=float(Tl), R=float(Rl), A=float(Al))
N['write_relax_time_critical_fs'] = hbar_eVs/(2*kappa)*1e15
N['mode_ratio_gamma'] = kappa/gA

# band-limited thermal injection from a 300 K contact (only |E|<2t propagates)
kT = kB*T0
E = np.linspace(-2*t+1e-6, 2*t-1e-6, 20001)
w = 1/(4*kT)/np.cosh(E/(2*kT))**2
w /= np.trapezoid(w, E)
_, _, Ath = scat(E, 0.0, kappa)
N['A_thermal_bandlimited_critical'] = float(np.trapezoid(Ath*w, E))
N['thermal_FWHM_meV'] = 3.5255*kT*1e3
N['band_width_meV'] = 4*t*1e3
N['resonance_FWHM_meV'] = (Gamma_el + 2*kappa)*1e3
# narrow injector (Gaussian, sigma) at resonance
for sig in (0.0005, 0.001, 0.002):
    wg = np.exp(-E**2/(2*sig**2)); wg /= np.trapezoid(wg, E)
    N[f'A_injector_sigma{sig*1e3:g}meV'] = float(np.trapezoid(Ath*wg, E))

# off cells: level pushed out of the band
for d in (0.3, 0.65, D_OFF):
    Tt, Rr, Aa = scat(0.0, d, gA)
    N[f'off_{d}eV'] = dict(T=float(Tt), R=float(Rr), A=float(Aa))
off = N[f'off_{D_OFF}eV']
N['tw_meV'] = tw*1e3
N['loss_per_off_cell'] = off['A'] + off['R']
N['wire_loss_127_off'] = 1 - (1 - N['loss_per_off_cell'])**(CELLS_PER_WIRE - 1)

# occupied cell: its next (1->2) level is E_C above; use EC = 0.25 eV (5-DB estimate)
EC_est = 0.25
Tt, Rr, Aa = scat(0.0, EC_est, gA)
N['occupied_EC0.25'] = dict(T=float(Tt), R=float(Rr), A=float(Aa))

# ------------------------------------------------------------------ 2. two-level cluster: capture level and storage level
def split(Ec):     # antibonding - bonding, centre energy Ec relative to outer sites
    return 2*math.sqrt((Ec/2)**2 + (2*tc)**2)
N['split_at_Ec0_meV'] = split(0)*1e3
Ec_res = 2*math.sqrt(max((hw_ph/2)**2 - (2*tc)**2, 0))
N['Ec_for_split_eq_phonon_meV'] = Ec_res*1e3
N['split_at_Ec0.15_meV'] = split(0.15)*1e3
# detailed balance: a stored electron in an addressed column escapes upward at
# gamma_up = gamma_down * exp(-Ds/kT). For write-mode relaxation tau_down = 0.1 ps,
# an exposure window W and loss <= 1e-9 per exposure: Ds >= kT ln(W/(1e-9 tau_down))
tau_down = 0.1e-12
for W in (100e-12, 1e-9):
    N[f'Ds_min_W{W*1e12:g}ps_eV'] = kB*T0*math.log(W/(1e-9*tau_down))
Ds = N['Ds_min_W100ps_eV']
N['Ec_for_Ds_min_eV'] = math.sqrt(max(Ds**2 - 4*(2*tc)**2, 0))
N['gate_for_Ds_min_V'] = N['Ec_for_Ds_min_eV']/alpha_g
N['escape_time_at_60meV_ps'] = tau_down*math.exp(0.060/(kB*T0))*1e12

# ------------------------------------------------------------------ 3. retention with BCH codes
def tau(Eb, nu=nu, T=T0):
    return math.exp(Eb/(kB*T))/nu

def min_bch(Eb, n_bits, dt=scrub, T=T0, nu=nu, tmax=20):
    p = dt/tau(Eb, nu, T)
    nw = n_bits/64
    for tt in range(1, tmax+1):
        n = 64 + 7*tt
        rate = nw*math.comb(n, tt+1)*p**(tt+1)/dt
        if rate*YEAR < 1:
            return tt, n
    return None, None

def Eb_req(n_bits, tt, dt=scrub, nu=nu, T=T0):
    n = 64 + 7*tt
    nw = n_bits/64
    # nw*C(n,t+1)*p^(t+1)/dt < 1/yr
    p = (dt/(YEAR*nw*math.comb(n, tt+1)))**(1/(tt+1))
    return kB*T*math.log(nu*dt/p)

# ------------------------------------------------------------------ 4. dual-rail crosstalk and pitch
def shift_monopole(px, py, eps=5.6, lam=5.0, R=60.0):
    nx, ny = int(R/px)+1, int(R/py)+1
    X, Y = np.meshgrid(np.arange(-nx, nx+1)*px, np.arange(-ny, ny+1)*py, indexing='ij')
    r = np.hypot(X, Y); m = (r > 0) & (r <= R)
    return float(np.sum(COUL/(eps*r[m])*np.exp(-r[m]/lam)))

def shift_dualrail(px, py, eps=5.6, lam=5.0, R=60.0):
    """Bits are pairs of cells (2m px, n py),(2m+1 px, n py). Mean occupation 1/2 is
    compensated; the target's own partner is deterministic. Worst case: every other
    pair oriented to add."""
    V = lambda r: COUL/(eps*r)*np.exp(-r/lam)
    nx, ny = int(R/(2*px))+1, int(R/py)+1
    tot = 0.0
    for m in range(-nx, nx+1):
        for n in range(-ny, ny+1):
            if m == 0 and n == 0:
                continue
            r1 = math.hypot(2*m*px, n*py); r2 = math.hypot((2*m+1)*px, n*py)
            if min(r1, r2) > R:
                continue
            tot += abs(V(r1) - V(r2))/2
    return tot

def min_pitch(fn, crit, pxmin=3, pymin=4, **kw):
    best = None
    for m in range(pxmin, 80):
        for n in range(pymin, 40):
            px, py = m*a, n*b
            if best is not None and px*py >= best[0]*best[1]:
                break
            if fn(px, py, **kw) <= crit:
                best = (px, py); break
    return best

crit = Gamma_el/2      # keep resonance within half its elastic width
N['crit_meV'] = crit*1e3
N['pitch_dualrail'] = {lam: min_pitch(shift_dualrail, crit, lam=lam) for lam in (2.0, 5.0, 10.0)}
N['pitch_monopole'] = {lam: min_pitch(shift_monopole, crit, lam=lam) for lam in (2.0, 5.0, 10.0)}
N['shift_old_pitch_monopole'] = shift_monopole(1.152, 1.152)

px, py = N['pitch_dualrail'][5.0]
cells_cm2 = 1/(px*py*1e-14)
bits = cells_cm2*area/2
N['design'] = dict(px=px, py=py, cells_cm2=cells_cm2, logical_bits=bits, TB=bits/8e12, GB=bits/8e9)
pxm, pym = N['pitch_monopole'][5.0]
N['design_monopole_GB'] = 1/(pxm*pym*1e-14)*area/8e9

# retention requirements for the design array
N['retention'] = {}
for tt in (1, 3, 5):
    N['retention'][tt] = dict(n=64+7*tt, overhead=7*tt/64,
                              Eb=[Eb_req(bits, tt, nu=x) for x in (nu_lo, nu_hi)])
N['retention_noECC'] = [kB*T0*math.log(x*bits*YEAR) for x in (nu_lo, nu_hi)]
N['bch_vs_Eb'] = {Eb: min_bch(Eb, bits) for Eb in (0.8, 0.9, 1.0, 1.05, 1.1)}
N['tau_0.65_ms'] = tau(0.65)*1e3
N['tau_1.0_s'] = tau(1.0)
N['tau_ratio_300_350_at1.0'] = tau(1.0)/tau(1.0, T=350)
N['bch_at_350K_Eb1.0'] = min_bch(1.0, bits, T=350)
N['Eb_req_350K_t5'] = Eb_req(bits, 5, T=350)
N['Tmax_Eb1.05_t5'] = next(T for T in np.arange(300, 400, 0.5) if Eb_req(bits, 5, T=T) > 1.05)
N['hold_gate_swing_V'] = 1.0/alpha_g
Si_gap = 1.12
N['Si_gap'] = Si_gap

# ------------------------------------------------------------------ 5. wire topology and throughput
vg = 2*t*e*a*1e-9/hbar
N['v_g'] = vg
cells_per_word = 128
word_len_nm = cells_per_word*px
N['word_length_nm'] = word_len_nm
words_per_wire = 1
N['words_per_wire'] = words_per_wire
L_nm = words_per_wire*word_len_nm
N['wire_length_um'] = L_nm/1e3
N['transit_ps'] = L_nm*1e-9/vg*1e12
rows_cm = 1/(py*1e-7)
n_wires = rows_cm*(math.sqrt(area)/(L_nm*1e-7))*math.sqrt(area)
N['wires'] = n_wires
N['injections_per_s'] = n_wires/(N['transit_ps']*1e-12)
E_inj = 1.0  # eV above contact Fermi level during operation (>= hold barrier)
N['energy_per_injection_J'] = E_inj*e
N['dataplane_power_all_wires_W'] = N['injections_per_s']*E_inj*e
# dual-rail NAND: 3 injections (one AND into m-, two ORs into m+), plus 2 resets
N['nand_injections'] = 3
N['energy_nand_J'] = 3*E_inj*e
# attempts at critical capture 0.5 with verify
q = 1 - N['critical']['A']
N['attempts_bit_9999'] = next(n for n in range(1, 99) if 1 - q**n >= 0.9999)

# ------------------------------------------------------------------ 5b. customised Blocks (FZC in the fabric)
# source Block: releases its stored electron into the wire; energy spread = its own
# coupling width to the wire, chosen <= 1 meV (Section 2.4)
G_src = 0.001
N['source_width_meV'] = G_src*1e3
N['source_release_time_fs'] = hbar_eVs/G_src*1e15
wg = np.exp(-E**2/(2*(G_src/2.355)**2)); wg /= np.trapezoid(wg, E)
# Lorentzian source line of FWHM G_src
wl = (G_src/2)**2/(E**2 + (G_src/2)**2); wl /= np.trapezoid(wl, E)
N['A_from_source_block'] = float(np.trapezoid(Ath*wl, E))
# control Block: an excess electron at distance r shifts a neighbour by V(r)
def Vscr(r, eps=5.6, lam=5.0):
    return COUL/(eps*r)*math.exp(-r/lam)
N['control_shift_eV'] = {r: Vscr(r) for r in (0.768, 1.152, 1.536)}
# predication: a source Block shifted by V no longer releases on the resonant pulse;
# residual release ~ Lorentzian tail of its own line
for r in (0.768, 1.152):
    d = Vscr(r)
    N[f'predicate_leak_r{r}'] = (G_src/2)**2/(d**2 + (G_src/2)**2)
# edge Blocks per wire: source + control at one end, sink at the other
N['edge_blocks_per_wire'] = 3
N['edge_overhead'] = 3/CELLS_PER_WIRE
# shift a control Block induces on the nearest data cell one tile pitch away
N['control_to_data_shift_meV'] = {k: Vscr(k*px)*1e3 for k in (1, 2)}

# ------------------------------------------------------------------ 6. fabrication and stability
atoms = cells_cm2*area*5 + n_wires*L_nm/a
N['atoms'] = atoms
N['fab_serial_years'] = atoms/1e4/YEAR
N['fab_areal_years'] = area*1e8/0.053/24/365.25
def Edes(frac, years, nu=1e13):
    return kB*T0*math.log(nu*years*YEAR/(-math.log(1 - frac)))
N['Edes_1e-6_10yr'] = Edes(1e-6, 10)

json.dump(N, open('numbers_r5.json', 'w'), indent=1, default=lambda o: float(o) if hasattr(o, '__float__') else str(o))

if '--figs' in sys.argv:
    import fea_r5_figs
    fea_r5_figs.make_all(N, globals())

if __name__ == '__main__':
    for k, v in N.items():
        print(k, '=', v)
