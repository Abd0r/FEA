#!/usr/bin/env python3
"""FEA R6: three-tier hybrid (electron working tier + hydrogen-position memory).
Reuses the R5 capture model (fea_r5.py). Writes numbers_r6.json. Use --figs for figures."""
import io, json, math, sys, contextlib, importlib.util
import numpy as np

spec = importlib.util.spec_from_file_location('r5', 'fea_r5.py')
r5 = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(r5)

kB, YEAR, COUL = r5.kB, r5.YEAR, r5.COUL
a, b, c = r5.a, r5.b, r5.c
kappa, scat = r5.kappa, r5.scat
N = {}

def V(r, eps=5.6, lam=5.0):
    return COUL/(eps*r)*math.exp(-r/lam)

# ---------------------------------------------------------------- 1. hydrogen-position memory tier
# measured (STM atom tracking, Si(001)): intradimer hop 1.01 +/- 0.05 eV, nu = 10^10.3 Hz;
# intrarow hop 1.75 +/- 0.10 eV, nu = 10^14.5 Hz
H = {'dimer': (1.01, 0.05, 10**10.3), 'row': (1.75, 0.10, 10**14.5)}
def tau(E, nu, T=300.0):
    return math.exp(E/(kB*T))/nu
for k, (E, dE, nu) in H.items():
    N[f'tau_{k}'] = {T: tau(E, nu, T) for T in (300, 330, 350)}
    N[f'tau_{k}_range_300'] = (tau(E - dE, nu), tau(E + dE, nu))

def bch_needed(tau_s, n_bits, dt):
    p = dt/tau_s; nw = n_bits/64
    for t in range(1, 30):
        n = 64 + (8 if t == 1 else 7*t)
        if nw*math.comb(n, t + 1)*p**(t + 1)/dt*YEAR < 1:
            return t, n, p
    return None, None, p

# ---------------------------------------------------------------- 2. working tier (electron) estimate
Eg2 = 0.56
N['working_tau_s'] = (tau(Eg2, 15.6e12), tau(Eg2, 1e13))

# ---------------------------------------------------------------- 3. reading a hydrogen bit through its port Block
# the DB of the memory dimer sits at distance r (H on far atom) or r - c (H on near atom)
def contrast(r):
    return V(r - c) - V(r)
for r in (0.768, 1.0, 1.2, 1.536):
    N[f'read_contrast_r{r}_meV'] = contrast(r)*1e3
gL = r5.N['logic_gammaA_over_kappa']*kappa
def read_error(r):
    # logic-mode port tuned to resonance for memory state 0 (reflects);
    # state 1 detunes it by the contrast; error = reflection of the detuned port
    Tt, Rr, Aa = scat(0.0, contrast(r), gL)
    return float(Rr + Aa)
r_port = 0.768
N['r_port_nm'] = r_port
N['read_err_single'] = read_error(r_port)
N['read_err_r1.2'] = read_error(1.2)
e1 = N['read_err_single']
N['read_err_majority3'] = 3*e1**2 - 2*e1**3
N['read_err_majority5'] = sum(math.comb(5, j)*e1**j*(1 - e1)**(5 - j) for j in (3, 4, 5))

# ---------------------------------------------------------------- 4. symmetric dual rail across the wire
# rail+ at (0, +d), rail- at (0, -d) with the wire on y = 0: the potential on every
# wire site is identical whichever rail is charged.
d_cell = 2*b
def wire_asym(dx_err):
    # residual data-dependent potential on the wire if one rail is displaced by dx_err along x
    return abs(V(math.hypot(dx_err, d_cell)) - V(d_cell))
N['wire_static_potential_meV'] = V(d_cell)*1e3
N['wire_asym_one_site_meV'] = wire_asym(a)*1e3

# ---------------------------------------------------------------- 5. crosstalk at a working Block, all sources
def shift_total(px, py, eps=5.6, lam=5.0, R=60.0):
    """Worst-case data-dependent shift at a working Block at (0, +d) from
    (i) other dual-rail working pairs, rails at (x, y0 +/- d), and
    (ii) other dual-rail hydrogen-memory pairs: a dimer beside each rail Block,
    mirrored across the wire and holding complementary H positions.
    The target's own pair and own memory pair are excluded (deterministic / signal)."""
    Vl = lambda r: V(r, eps, lam)
    tot = 0.0
    nx, ny = int(R/px) + 1, int(R/py) + 1
    tx, ty = 0.0, d_cell
    dist = lambda x, y: math.hypot(x - tx, y - ty)
    D, Dn = d_cell + r_port, d_cell + r_port - c       # DB far / near positions
    for i in range(-nx, nx + 1):
        for j in range(-ny, ny + 1):
            if i == 0 and j == 0:
                continue
            x, y0 = i*px, j*py
            r1, r2 = dist(x, y0 + d_cell), dist(x, y0 - d_cell)
            if min(r1, r2) > R:
                continue
            tot += abs(Vl(r1) - Vl(r2))/2
            sA = Vl(dist(x, y0 + Dn)) + Vl(dist(x, y0 - D))
            sB = Vl(dist(x, y0 + D)) + Vl(dist(x, y0 - Dn))
            tot += abs(sA - sB)/2
    return tot

def wire_memory_residual():
    # wire site at origin: mirrored complementary memory pair -> identical in both states
    D, Dn = d_cell + r_port, d_cell + r_port - c
    return abs((V(Dn) + V(D)) - (V(D) + V(Dn)))
N['wire_memory_residual_meV'] = wire_memory_residual()*1e3
N['wire_single_memory_bit_meV'] = (V(d_cell + r_port - c) - V(d_cell + r_port))*1e3

crit = kappa
best = None
for m in range(3, 60):
    for n in range(6, 40):
        px, py = m*a, n*b
        if best is not None and px*py >= best[0]*best[1]:
            break
        if shift_total(px, py) <= crit:
            best = (px, py); break
N['pitch'] = best
px, py = best
pairs_cm2 = 1/(px*py*1e-14)
area = 0.5
N['work_bits'] = pairs_cm2*area
N['mem_bits'] = pairs_cm2*area             # one dual-rail hydrogen bit per Block pair
N['mem_GB'] = N['mem_bits']/8e9
N['work_GB'] = N['work_bits']/8e9
for lam in (2.0, 10.0):
    bb = None
    for m in range(3, 80):
        for n in range(6, 50):
            qx, qy = m*a, n*b
            if bb is not None and qx*qy >= bb[0]*bb[1]:
                break
            if shift_total(qx, qy, lam=lam) <= crit:
                bb = (qx, qy); break
    N[f'pitch_lam{lam}'] = bb
    N[f'mem_GB_lam{lam}'] = 1/(bb[0]*bb[1]*1e-14)*area/8e9

# BCH for the dimer tier at the computed capacity
for dt in (60, 3600):
    t, n, p = bch_needed(N['tau_dimer'][300], N['mem_bits'], dt)
    N[f'dimer_bch_dt{dt}'] = dict(t=t, n=n, p=p, overhead=(n - 64)/64 if n else None)
t, n, p = bch_needed(N['tau_dimer'][330], N['mem_bits'], 60)
N['dimer_bch_330K_dt60'] = dict(t=t, n=n, p=p)
t, n, p = bch_needed(tau(1.01 - 0.05, 10**10.3), N['mem_bits'], 60)
N['dimer_bch_low_barrier_dt60'] = dict(t=t, n=n, p=p)

# scrub load: read every memory bit each minute through ports (majority of 3 reads) ->
# injections per wire per second, vs capacity
word_cells = 128
wire_len = 64*px              # 64 pairs along the wire
N['wire_len_um'] = wire_len/1e3
transit = wire_len*1e-9/r5.vg
N['transit_ps'] = transit*1e12
mem_bits_per_wire = 128
inj_per_s_scrub = mem_bits_per_wire*3/60
N['scrub_fraction'] = inj_per_s_scrub*transit

# ---------------------------------------------------------------- 5b. compensation, wires, throughput, lines
# static (data-independent) modulation of the wire potential by the 2 electrons of each pair
mod = 2*V(d_cell) - 2*V(math.hypot(px/2, d_cell))
N['wire_modulation_meV'] = mod*1e3
N['compensation_charge_depth_nm'] = d_cell
pairs = pairs_cm2*area
wires = pairs/64
N['wires'] = wires
N['injections_per_s'] = wires/transit
N['power_full_W'] = N['injections_per_s']*1.602e-19
# NAND: three write steps, each repeated until verified (expected 1/A attempts + one verify read)
A = 0.5
N['nand_expected_injections'] = 3*(1/A + 1)
N['nand_energy_J'] = N['nand_expected_injections']*1.602e-19
N['operand_disturb_per_nand'] = N['nand_expected_injections']*1e-6
# column lines: assumed sheet resistance (not measured for this geometry)
Rs, wl, Lsub, Cpl, Vsw = 1e3, 50e-9, 100e-6, 0.2e-15/1e-6, 3.3
R = Rs*Lsub/wl; C = Cpl*Lsub
N['line'] = dict(R_ohm=R, C_F=C, RC_s=R*C, E_switch_J=C*Vsw**2, wires_per_line=Lsub/(py*1e-9),
                 E_per_wire_op_J=C*Vsw**2/(Lsub/(py*1e-9)))
# control Block two rows from the source Block
N['control_shift_1.536_eV'] = V(1.536)
N['predicate_leak_1.536'] = (0.0005)**2/(V(1.536)**2 + 0.0005**2)
# patterned atoms: 2 Blocks x5 + 2 memory DBs per pair + wire (2 atoms per a, bare-dimer line)
N['atoms'] = pairs*(10 + 2 + 2*px/a)
N['fab_serial_years'] = N['atoms']/1e4/YEAR

# ---------------------------------------------------------------- 6. write cost of a hydrogen bit (unknown yield)
for Y in (1e-4, 1e-6, 1e-8):
    N[f'write_energy_J_yield{Y:g}'] = 2.0*1.602e-19/Y   # ~2 eV electrons
N['work_energy_inj_J'] = 1.602e-19

json.dump(N, open('numbers_r6.json', 'w'), indent=1, default=lambda o: float(o) if hasattr(o, '__float__') else str(o))
if '--figs' in sys.argv:
    import fea_r6_figs
    fea_r6_figs.make_all(N, globals())
if __name__ == '__main__':
    for k, v in N.items():
        print(k, '=', v)
