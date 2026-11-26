#!/usr/bin/env python3
"""FEA R7: room-temperature compute + working memory (no hydrogen storage tier).
Every number in FEA-R7.tex comes from this script -> numbers_r7.json.

Parts
  A  chip: crosstalk-limited pitch WITHOUT hydrogen memory dimers, capacity, read model, energies
  B  DFT (GPAW, ../dft/runs): charge-resolved H-hop barriers = atomic-pattern stability, DB- mechanism,
     neighbour-charge tilt vs the point-charge model, charge relaxation (lambda) when cc.json exists
  C  retention + chip simulation: imported from ../opt (they read the chip block written here)
Reuses the R5 capture model (fea_r5.py).
"""
import io, json, math, os, sys, contextlib, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
spec = importlib.util.spec_from_file_location('r5', 'fea_r5.py')
r5 = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(r5)

kB, YEAR, COUL = r5.kB, r5.YEAR, r5.COUL
a, b, c = r5.a, r5.b, r5.c
kappa = r5.kappa
N = {}
AREA = 0.5                     # cm^2 reference die
d_cell = 2*b                   # rails at +-1.536 nm from the wire

def V(r, eps=5.6, lam=5.0):    # screened Coulomb, eV, r in nm
    return COUL/(eps*r)*math.exp(-r/lam)

# ------------------------------------------------------------------ A1 crosstalk without memory dimers
def shift_total(px, py, eps=5.6, lam=5.0, R=60.0):
    """Worst-case data-dependent shift at a working Block at (0,+d) from all other dual-rail pairs."""
    tot = 0.0; nx, ny = int(R/px) + 1, int(R/py) + 1
    for i in range(-nx, nx + 1):
        for j in range(-ny, ny + 1):
            if i == 0 and j == 0: continue
            x, y0 = i*px, j*py
            r1, r2 = math.hypot(x, y0 + d_cell - d_cell), math.hypot(x, y0 - d_cell - d_cell)
            if min(r1, r2) > R: continue
            tot += abs(V(r1, eps, lam) - V(r2, eps, lam))/2
    return tot

def best_pitch(lam):
    best = None
    for m in range(3, 80):
        for n in range(6, 50):
            px, py = m*a, n*b
            if best is not None and px*py >= best[0]*best[1]: break
            if shift_total(px, py, lam=lam) <= kappa:
                best = (px, py); break
    return best

N['crit_meV'] = kappa*1e3
for lam in (2.0, 5.0, 10.0):
    p = best_pitch(lam)
    bits = AREA/(p[0]*p[1]*1e-14)
    N[f'pitch_lam{lam}'] = p
    N[f'raw_GB_lam{lam}'] = bits/8e9
P = N['pitch_lam5.0']
RAW_BITS = AREA/(P[0]*P[1]*1e-14)

# ------------------------------------------------------------------ A2 read of an electron bit (R5 logic mode)
LM = r5.N['logic_mode']                    # at the D3 operating point (gamma_A/kappa = 5e-7, 0.23 us)
N['read_silent_err'] = LM['T']             # empty Block passes -> misread (silent)
N['read_disturb'] = LM['A']                # empty Block captures during read -> '11' in dual rail = erasure

# ------------------------------------------------------------------ A3 energies (R6 line model, R5 injection)
INJ = 1.602176634e-19
LINE_J = 0.2178e-12/6200.396825396825      # R6: 0.22 pJ per line switch shared by 6,200 wires -> J per wire op
N['E_injection_aJ'] = INJ*1e18
# A buried column line sets the mode of ONE bit position across the 6,200 Words it crosses, so its switching
# energy is paid per bit operated (per wire per bit-step), not shared over the bits of a Word.
# Line requirement L1 (optimisation C): 10x cheaper switching than the R6 line (3.3 V swing), i.e. a swing of
# 3.3/sqrt(10) V. Detuning a Block by the storage depth Ds then needs a lever arm Ds/swing.
R5N = r5.N
V_R6 = 3.3; LINE_SCALE = 0.1
# Optimisation E2: the edge drivers charge each column line in LINE_STEPS equal voltage steps from tank
# capacitors (stepwise/adiabatic charging) and recover the charge on the way down, so the energy dissipated
# per switch is CV^2/(2*LINE_STEPS) instead of CV^2/2. Each step needs a few RC, so a switch takes ~LINE_STEPS x longer.
LINE_STEPS = 4
N['line'] = dict(E_R6_aJ=LINE_J*1e18, E_req_aJ=LINE_J*LINE_SCALE*1e18, swing_R6_V=V_R6,
                 swing_req_V=V_R6*math.sqrt(LINE_SCALE), Ds_eV=R5N['Ds_min_W100ps_eV'])
N['line']['lever_req_eV_per_V'] = N['line']['Ds_eV']/N['line']['swing_req_V']
N['line']['lever_R5_eV_per_V'] = 0.30
N['line'].update(steps=LINE_STEPS, E_stepped_aJ=LINE_J*LINE_SCALE/LINE_STEPS*1e18,
                  swing_same_energy_V=V_R6*math.sqrt(LINE_SCALE*LINE_STEPS))   # alternative use: keep 3.5 aJ, relax the swing
N['line']['lever_same_energy_eV_per_V'] = N['line']['Ds_eV']/N['line']['swing_same_energy_V']
def energies(line_J):
    # refresh / access of one bit: one logic-mode read (1 injection; rewrites are rare) + 2 line switches
    # (column into logic mode and back to hold); bitwise NAND: 9 injections + 3 write steps x 2 line switches
    return dict(refresh=INJ + 2*line_J, gate=9*INJ + 6*line_J)
E_req, E_r6 = energies(LINE_J*LINE_SCALE/LINE_STEPS), energies(LINE_J)
E_conv = energies(LINE_J*LINE_SCALE)                                   # L1 lines, conventional (one-step) drive
N['E_line_per_wire_op_aJ'] = LINE_J*LINE_SCALE/LINE_STEPS*1e18
N['E_refresh_bit_conv_aJ'] = E_conv['refresh']*1e18; N['E_gate_conv_aJ'] = E_conv['gate']*1e18
N['E_refresh_bit_aJ'] = E_req['refresh']*1e18; N['E_gate_aJ'] = E_req['gate']*1e18
N['E_refresh_bit_R6lines_aJ'] = E_r6['refresh']*1e18; N['E_gate_R6lines_aJ'] = E_r6['gate']*1e18
N['E_access_bit_aJ'] = N['E_refresh_bit_aJ']

# ------------------------------------------------------------------ A4 code and chip block for the simulators
# Optimisation B: a shortened binary BCH code of length 64 + 7t with t = 8 (minimum distance 17). Dual-rail
# turns lost electrons into erasures, so it corrects 16 erasures (or 8 errors) per Word.
CODE_T = 8; CODE_N = 64 + 7*CODE_T; CODE_D = 2*CODE_T + 1
# Optimisation E1: the FZC decodes erasures by Gaussian elimination (maximum likelihood). Any s <= d-1 = 16
# erasures decode; beyond that a pattern fails when a codeword hides inside it. Monte Carlo on the real
# shortened BCH(120,64) (opt/bch_erasure.py) gives no failure up to s = 44 and agrees with the random-code
# value 2^-(n-k-s) (or is below it) from s = 45 to 56; the simulators use 2^-(n-k-s).
N['code'] = dict(t=CODE_T, n=CODE_N, d=CODE_D, k=64, overhead=(CODE_N - 64)/64, decoder='ML')
_be = os.path.join(HERE, '..', 'opt', 'bch_erasure.json')
if os.path.exists(_be):
    _rows = json.load(open(_be))['rows']
    N['code']['ml_check'] = dict(trials=_rows[-1]['trials'],
        s_max_no_fail=max(r['s'] for r in _rows if r['n_fail'] == 0 and all(q['n_fail'] == 0 for q in _rows if q['s'] <= r['s'])),
        p_fail_s48=next(r['p_fail'] for r in _rows if r['s'] == 48), rc_s48=2.0**-(CODE_N - 64 - 48))
words = RAW_BITS/CODE_N
N['chip'] = dict(area_cm2=AREA, pitch_nm=P, raw_bits=RAW_BITS, raw_GB=RAW_BITS/8e9, words=words,
                 usable_GB=words*64/8e9, read_silent_err=N['read_silent_err'], read_disturb=N['read_disturb'],
                 code_n=CODE_N, code_d=CODE_D, code_k=64, decoder='ML',
                 E_refresh_J=E_req['refresh'], E_gate_J=E_req['gate'], E_access_J=E_req['refresh'],
                 E_refresh_R6lines_J=E_r6['refresh'], E_gate_R6lines_J=E_r6['gate'], line_RC_s=40e-9,
                 E_refresh_conv_J=E_conv['refresh'], line_steps=LINE_STEPS)
N['usable_GB_code92'] = RAW_BITS/92*64/8e9


# ------------------------------------------------------------------ B  DFT results (GPAW PBE-DZP, ../dft/runs)
RUNS = os.path.join(HERE, '..', 'dft', 'runs')
def J(name):
    f = os.path.join(RUNS, name); return json.load(open(f)) if os.path.exists(f) else None

def band_geometry(fname, target=(89, 93)):
    """Per image: Si_i-H, Si_k-H, dimer length, energy relative to image 0."""
    from ase.io import read
    imgs = read(os.path.join(RUNS, fname), index=':'); i, k = target; out = []
    E0 = imgs[0].get_potential_energy()
    for im in imgs:
        h = [n for n, x in enumerate(im.get_chemical_symbols()) if x == 'H']
        di = im.get_distances(i, h, mic=True); j = h[int(np.argmin(di))]
        out.append(dict(SiH_i=float(di.min()), SiH_k=float(im.get_distance(k, j, mic=True)),
                        dimer=float(im.get_distance(i, k, mic=True)), dE=float(im.get_potential_energy() - E0)))
    return out

NEB = {'DB0': 'neb_q0', 'DB-': 'neb_qm5', 'DB+': 'neb_qp'}
N['dft'] = {}
for lab, f in NEB.items():
    r = J(f + '.json')
    N['dft'][lab] = dict(Ea=r['Ea_fwd'], steps=r['steps'], converged=r['converged'], n_images=len(r['energies']) - 2,
                         profile=band_geometry(f + '_band.extxyz'))
prof = N['dft']['DB-']['profile']
mid = len(prof)//2
N['dft']['DB-']['bridge_dE'] = prof[mid]['dE']
N['dft']['DB-']['bridge_dimer'] = prof[mid]['dimer']
N['dft']['DB-']['dimer_endpoint'] = prof[0]['dimer']
sp = J('spin_fix.json')
N['dft']['spin_check'] = dict(Ea_spin=sp['Ea_spin_corrected'], magmom_min=sp['points']['0']['magmom'],
                              magmom_saddle=sp['points'][str(sp['saddle_image'])]['magmom'])

# lattice relaxation on charging: displacement of the DB atom (index 93 carries the DB in the L structures)
from ase.io import read as _read
_g = {q: _read(os.path.join(RUNS, f'{q}_L.extxyz')) for q in ('q0', 'qm', 'qp')}
def _dz(q):
    d = _g[q].positions[93] - _g['q0'].positions[93]; return float(d[2])
N['dft']['DB_atom_dz_on_charging'] = {'DB-': _dz('qm'), 'DB+': _dz('qp')}

# pattern stability: rate of one H hop per DB site (Arrhenius), prefactor range
NU = {'Hill_clean_meas': 10**10.3, 'Hill_low': 10**9.8, 'Hill_high': 10**10.8, 'phonon': 1e13}
def hop_rate(Ea, nu, T=300.0, two_saddles=False):
    r = nu*math.exp(-Ea/(kB*T)); return r/2 if two_saddles else r
N['stability'] = {}
for lab in ('DB0', 'DB-', 'DB+'):
    Ea = N['dft'][lab]['Ea']; two = lab == 'DB-'   # DB- passes a bridged intermediate: half its arrivals return
    N['stability'][lab] = {k: dict(rate_per_s=hop_rate(Ea, nu, two_saddles=two),
                                   per_DB_year=hop_rate(Ea, nu, two_saddles=two)*YEAR,
                                   mean_life_yr=1/hop_rate(Ea, nu, two_saddles=two)/YEAR) for k, nu in NU.items()}
DB_PER_BIT = 10                                   # two five-DB Blocks per dual-rail bit
N['n_DB_chip'] = RAW_BITS*DB_PER_BIT
N['pattern_defects_per_year_chip'] = {k: N['stability']['DB0'][k]['per_DB_year']*N['n_DB_chip'] for k in NU}
# DB+ design rule: positive occupancy fraction f+ at which DB+ hops equal DB0 hops
N['f_plus_max'] = {k: N['stability']['DB0'][k]['rate_per_s']/N['stability']['DB+'][k]['rate_per_s'] for k in NU}

# bonding / anneal limit: highest temperature for a 1 h anneal moving at most a fraction F of DBs
def T_max(Ea, nu, t, F):
    return Ea/(kB*math.log(nu*t/F))
N['anneal_Tmax_C'] = {f'{F:g}': {k: T_max(N['dft']['DB0']['Ea'], nu, 3600.0, F) - 273.15 for k, nu in NU.items()}
                      for F in (1e-6, 1e-4, 1e-2)}

# neighbour-charge tilt: DFT vs screened point charges at the relaxed distances
t2 = J('tilt_qm2_L.json'), J('tilt_qm2_R.json')
if all(t2):
    gL, gR = _read(os.path.join(RUNS, 'tilt_qm2_L.extxyz')), _read(os.path.join(RUNS, 'tilt_qm2_R.extxyz'))
    rn = float(gL.get_distance(93, 92, mic=True))/10; rf = float(gR.get_distance(89, 92, mic=True))/10   # nm
    N['tilt'] = dict(dE_DFT=t2[0]['energy_eV'] - t2[1]['energy_eV'], r_near_nm=rn, r_far_nm=rf,
                     dE_model_eps5_6=V(rn) - V(rf), dE_model_eps6_35=V(rn, 6.35) - V(rf, 6.35),
                     favoured='DB far from the neighbour')
    N['tilt']['bias_ratio_300K'] = math.exp(N['tilt']['dE_DFT']/(kB*300))
t0 = J('tilt_q0_L.json'), J('tilt_q0_R.json')
st = J('spin_tilt.json')
if all(t0):
    eL, eR = t0[0]['energy_eV'], t0[1]['energy_eV']
    if st:   # take each configuration in its lowest spin state (relaxed state vs the other spin state)
        N['tilt']['spin'] = dict(near_singlet=eL, near_triplet=st['L']['E'], far_triplet=eR, far_broken_sym=st['R']['E'],
                                 near_ST_gap_meV=(st['L']['E'] - eL)*1e3, far_ST_gap_meV=(st['R']['E'] - eR)*1e3)
        eL, eR = min(eL, st['L']['E']), min(eR, st['R']['E'])
    N['tilt']['dE_DFT_neutral'] = eL - eR
    # the data-dependent (charge-dependent) part is what the point-charge crosstalk model represents
    N['tilt']['dE_DFT_charge_part'] = N['tilt']['dE_DFT'] - N['tilt']['dE_DFT_neutral']
    N['tilt']['model_over_DFT_eps5_6'] = N['tilt']['dE_model_eps5_6']/N['tilt']['dE_DFT_charge_part']
    N['tilt']['model_over_DFT_eps6_35'] = N['tilt']['dE_model_eps6_35']/N['tilt']['dE_DFT_charge_part']

# charge relaxation (configuration coordinate) when cc.json is available
cc = J('cc.json')
if cc:
    N['cc'] = cc
    lam_m = 0.5*(cc['lam']['neutral, at DB- geometry'] + cc['lam']['DB- at neutral geometry'])
    lam_p = 0.5*(cc['lam']['neutral, at DB+ geometry'] + cc['lam']['DB+ at neutral geometry'])
    N['lambda_eV'] = dict(minus=lam_m, plus=lam_p)

if __name__ == '__main__':
    json.dump(N, open('numbers_r7.json', 'w'), indent=1, default=float)
    for k in ('pitch_lam2.0', 'pitch_lam5.0', 'pitch_lam10.0', 'raw_GB_lam2.0', 'raw_GB_lam5.0', 'raw_GB_lam10.0',
              'read_silent_err', 'read_disturb', 'E_refresh_bit_aJ', 'E_gate_aJ'):
        print(k, N[k])
    print('usable GB', N['chip']['usable_GB'])
