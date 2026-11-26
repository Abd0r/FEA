"""Harmonic transition-state prefactor for the intra-dimer H hop (Vineyard):
   nu = prod(freqs at minimum) / prod(real freqs at saddle), over the H atom and its two dimer Si atoms
   (partial Hessian, 9 modes; the saddle has one imaginary mode). Finite differences, delta = 0.01 A, 2-point.
Usage: gpaw python htst.py LABEL MIN.extxyz BAND.extxyz SADDLE_INDEX CHARGE [DELTA]
DELTA (A, default 0.01). For DELTA != 0.01 the outputs carry a suffix _dNN. Forces converged to 1e-4 eV/A and
the jellium slab is pinned to the undisplaced geometry, so the finite differences see no artificial force noise.
Writes runs/htst_LABEL.json. Each displacement is cached by ASE (runs/vib_LABEL_*), so a killed run resumes."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read
from ase.vibrations import Vibrations
from gpaw.mpi import world
from relax import make_calc
label, fmin, fband, isad, q = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
delta = float(sys.argv[6]) if len(sys.argv) > 6 else 0.01
if delta != 0.01:
    label = f'{label}_d{int(round(delta*100)):02d}'
tight = {} if delta == 0.01 else dict(convergence={'density': 1e-6, 'forces': 1e-4})

def subset(at):
    i, k = 89, 93
    h = [n for n, s in enumerate(at.get_chemical_symbols()) if s == 'H']
    dmin = [min(at.get_distance(i, j, mic=True), at.get_distance(k, j, mic=True)) for j in h]
    return [i, k, h[int(np.argmin(dmin))]]

out = {}
for tag, at in (('min', read(fmin)), ('saddle', read(fband, index=isad))):
    at = at.copy()
    if q == 0:   # neutral DB: start spin-polarized so it does not fall into the non-magnetic state
        m = np.zeros(len(at)); m[89] = m[93] = 0.5; at.set_initial_magnetic_moments(m)
    zj = (at.positions[:, 2].max() + 5.0, at.cell[2, 2] - 2.0)   # fixed: displacements must not move the jellium
    at.calc = make_calc(at, 'dzp', q, txt=f'runs/htst_{label}_{tag}.txt', zjell=zj, **tight)
    idx = subset(at)
    vib = Vibrations(at, indices=idx, name=f'runs/vib_{label}_{tag}', delta=delta, nfree=2)
    vib.run()
    f = vib.get_frequencies()                      # cm^-1, imaginary returned as complex
    out[tag] = dict(indices=idx, freqs_cm=[[float(np.real(x)), float(np.imag(x))] for x in f])
c = 2.99792458e10                                  # cm/s
fm = np.array([x[0] for x in out['min']['freqs_cm']])
fs = np.array([x[0] for x in out['saddle']['freqs_cm'] if x[1] == 0 and x[0] > 0])
n_imag = sum(1 for x in out['saddle']['freqs_cm'] if x[1] != 0)
if len(fs) == len(fm) - 1 and n_imag == 1:
    out['nu_Hz'] = float(c*np.prod(fm)/np.prod(fs))
else:
    out['nu_Hz'] = None
out['n_imag_saddle'] = n_imag
if world.rank == 0:
    json.dump(out, open(f'runs/htst_{label}.json', 'w'), indent=1)
    print(label, 'prefactor', out['nu_Hz'], 'imaginary modes at saddle', n_imag)
