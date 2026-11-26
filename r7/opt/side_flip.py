#!/usr/bin/env python3
"""Side-flip tolerance of FEA Blocks.

An intra-dimer H hop moves a DB to the other atom of its dimer (0.225 nm across the row).
With the HTST prefactor this happens about 10 times per DB-year at 82 C, so the design must
tolerate any combination of DB sides. This script asks whether it does, using the same screened
point-charge model that sets the array pitch (eps = 5.6, lambda = 5 nm, confirmed by DFT).

Geometry (nm): dimer rows along x (a = 0.384), rows b = 0.768 apart in y, wire = bare-dimer row at y = 0.
A bit = two five-DB Blocks mirrored across the wire, centres at y = +-2b (+-1.54 nm); Block sites:
centre (0, Y), along-row (+-a, Y), cross-row (0, Y +- b). Each DB sits on one atom of its dimer: y offset +-0.1125.
Bits tile at the R7 pitch (px along the wire, py across).

Tests (each Monte Carlo over random sides of every DB, data held fixed):
 1. absorber detuning: shift of the target absorber's level relative to all-nominal sides,
    against the resonance budget Gamma_el/2 = 2.8 meV that sets the pitch;
 2. wire disorder: transmission loss of a 64-bit Word's wire (tight binding, t = 20 meV) when the
    static modulation is compensated for nominal sides and the sides then randomise;
 3. the same with bare-dimer Block sites (charge centroid at the dimer centre: no side to flip).
Usage: python3 side_flip.py   -> side_flip.json
"""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COUL, EPS, LAM = 1.4399645, 5.6, 5.0
a, b, half = 0.384, 0.768, 0.1125
R7 = json.load(open(os.path.join(HERE, '..', 'paper', 'numbers_r7.json')))
R5 = json.load(open(os.path.join(HERE, '..', 'paper', 'numbers_r5.json')))
PX, PY = 7.68, 16.128
CRIT = R5['Gamma_el_meV']/2/1e3          # eV
rng = np.random.default_rng(7)

def V(r):
    return COUL/(EPS*r)*np.exp(-r/LAM)

GEOMS = {   # Block sites relative to its centre; y is measured away from the wire
    'orig': [[0, 0], [a, 0], [-a, 0], [0, b], [0, -b]],          # R7: one cross-row DB 0.77 nm from the wire
    'far':  [[0, 0], [a, 0], [-a, 0], [0, b], [a, b]],           # both cross-row DBs on the far side
    'line': [[0, 0], [a, 0], [-a, 0], [2*a, 0], [-2*a, 0]],      # five DBs along one row
}
P = dict(geom='orig', off=1.0, t=0.020)   # off: side offset as a fraction of 0.1125 nm (1 = single DB; bare dimer: buckling-set)

def block_dbs(cx, cy, sides):
    """positions of the five DBs of a Block centred at (cx, cy); sides in {+1,-1}"""
    S = np.array(GEOMS[P['geom']], dtype=float)
    p = np.empty_like(S); p[:, 0] = cx + S[:, 0]; p[:, 1] = cy + np.sign(cy)*S[:, 1]
    p[:, 1] += sides*half*P['off']
    return p

def charge_weights(mode):
    return {'centre': np.array([1, 0, 0, 0, 0.]), 'outer': np.array([0, .25, .25, .25, .25]),
            'even': np.full(5, 0.2)}[mode]

def layout(nx=7, ny=4):
    """bit centres around the target bit at the origin"""
    return [(i*PX, j*PY) for i in range(-nx, nx+1) for j in range(-ny, ny+1)]

def absorber_test(mode, n_mc=2000):
    bits = layout()
    Y = 2*b
    w = charge_weights(mode)
    data = rng.integers(0, 2, len(bits))
    t0 = bits.index((0.0, 0.0)); data[t0] = 0   # target bit: upper rail empty (being written), lower rail full
    def energy(sides_all, abs_side):
        ap = np.array([0.0, Y + abs_side*half*P['off']])
        E = 0.0
        for k, (cx, cy) in enumerate(bits):
            yc = cy - Y if data[k] == 0 else cy + Y
            p = block_dbs(cx, yc, sides_all[k])
            E += float(np.sum(w*V(np.hypot(*(p - ap).T))))
        return E
    nom = energy(np.ones((len(bits), 5)), 1)
    d = np.abs([energy(rng.choice([-1, 1], size=(len(bits), 5)), rng.choice([-1, 1])) - nom for _ in range(n_mc)])
    return dict(mean_meV=float(d.mean()*1e3), p99_meV=float(np.percentile(d, 99)*1e3),
                frac_over_budget=float(np.mean(d > CRIT)), budget_meV=CRIT*1e3)

def transmission(dU, t, E=0.0):
    k = math.acos(-E/(2*t))
    M = np.eye(2, dtype=complex)
    for u in dU:
        M = np.array([[-(E - u)/t, -1], [1, 0]], dtype=complex) @ M
    e = np.exp(1j*k)
    vin = np.array([1, e**-1]); vr = np.array([1, e])        # (psi_0, psi_-1): incident, reflected
    N = len(dU); out_t = np.array([e**N, e**(N-1)])
    sol = np.linalg.solve(np.column_stack([M @ vr, -out_t]), -(M @ vin))
    return float(abs(sol[1])**2)

def wire_test(mode, n_mc=100):
    """64-bit Word along x; wire sites every a; static potential from the 64 occupied Blocks,
    compensated for nominal sides, then sides randomised. Loss = worst over five injection energies."""
    nbits = 64
    xs = np.arange(-20*a, nbits*PX + 20*a, a)
    w = charge_weights(mode)
    Y = 2*b; t = P['t']
    data = rng.integers(0, 2, nbits)
    def pot(sides):
        U = np.zeros_like(xs)
        for k in range(nbits):
            p = block_dbs(k*PX + PX/2, -Y if data[k] == 0 else Y, sides[k])
            for j in range(5):
                U += w[j]*V(np.hypot(xs - p[j, 0], p[j, 1]))
        return U
    U0 = pot(np.ones((nbits, 5)))
    assert abs(transmission(np.zeros_like(xs), t) - 1) < 1e-9
    Es = np.linspace(-0.6, 0.6, 5)*2*t
    losses, rms = [], []
    for _ in range(n_mc):
        dU = pot(rng.choice([-1, 1], size=(nbits, 5))) - U0
        rms.append(np.std(dU))
        losses.append(max(1 - transmission(dU, t, E) for E in Es))
    losses = np.array(losses)
    return dict(rms_dU_meV=float(np.mean(rms)*1e3), loss_median=float(np.median(losses)),
                loss_p99=float(np.percentile(losses, 99)), pass_frac=float(np.mean(losses <= 1e-3)))

if __name__ == '__main__':
    out = {'geometry': dict(px=PX, py=PY, rail_offset_nm=2*b, flip_nm=2*half, budget_meV=CRIT*1e3)}
    # 1. R7 Block, single DBs, t = 20 meV (the failing baseline)
    out['baseline'] = dict(absorber_centre=absorber_test('centre'), absorber_outer=absorber_test('outer'),
                           wire_centre=wire_test('centre'), wire_outer=wire_test('outer'))
    print('baseline', out['baseline'], flush=True)
    # 2. sweep: geometry x wire hopping x side offset (1 = single DB; <1 = residual buckling offset of a bare dimer)
    out['sweep'] = []
    for g in ('orig', 'far', 'line'):
        for off in (1.0, 0.3, 0.1):
            P.update(geom=g, off=off)
            ab = absorber_test('even', n_mc=400)
            for t in (0.020, 0.050, 0.100, 0.200):
                P['t'] = t
                wt = wire_test('even', n_mc=60)
                row = dict(geom=g, off=off, t_meV=t*1e3, absorber_p99_meV=ab['p99_meV'],
                           absorber_over=ab['frac_over_budget'], **wt)
                out['sweep'].append(row); print(row, flush=True)
    P.update(geom='orig', off=1.0, t=0.020)
    r = 4*b
    out['single_flip_partner_meV'] = float((V(r - 2*half) - V(r))*1e3)
    json.dump(out, open(os.path.join(HERE, 'side_flip.json'), 'w'), indent=1)
