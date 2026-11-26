"""Bare (fully depassivated) dimer: is its charge centroid on the dimer centre, or does it buckle?
A side-tolerant Block uses bare dimers as sites: with no H left in the dimer there is no intra-dimer hop,
and if the dimer stays symmetric (or flips much faster than the electron transit) its charge sits at the centre.
Usage: gpaw python bare.py START(sym|buck) CHARGE(0|-1)
Writes runs/bare_{start}_{q0|qm}.json with energy, buckling dz = z(89)-z(93), dimer length."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read, write
from ase.optimize import BFGS
from gpaw.mpi import world
from relax import make_calc
start, q = sys.argv[1], float(sys.argv[2])
tag = f"bare_{start}_{'q0' if q == 0 else 'qm'}"
at = read('slab_4x4x6_L.extxyz')
i, k = at.info['target_dimer']
hs = [n for n, s in enumerate(at.get_chemical_symbols()) if s == 'H']
hd = [n for n in hs if min(at.get_distance(i, n, mic=True), at.get_distance(k, n, mic=True)) < 1.8]
del at[hd]                                   # H atoms follow the Si atoms, so i, k are unchanged
assert at[i].symbol == 'Si' and at[k].symbol == 'Si'
if start == 'buck':
    at.positions[i, 2] += 0.35; at.positions[k, 2] -= 0.35
m = np.zeros(len(at))
if q != 0:
    m[i] = m[k] = 0.5
at.set_initial_magnetic_moments(m)
zj = (at.positions[:, 2].max() + 5.0, at.cell[2, 2] - 2.0)
at.calc = make_calc(at, 'dzp', q, txt=f'runs/{tag}.txt', zjell=zj)
opt = BFGS(at, trajectory=f'runs/{tag}.traj', logfile=f'runs/{tag}.log')
opt.run(fmax=0.03, steps=300)
e = at.get_potential_energy(); f = at.get_forces()
mag = float(at.calc.get_magnetic_moment()) if q != 0 else 0.0
write(f'runs/{tag}.extxyz', at)
if world.rank == 0:
    json.dump(dict(energy_eV=e, fmax=float(np.linalg.norm(f, axis=1).max()), steps=opt.nsteps,
                   dz_buckle=float(at.positions[i, 2] - at.positions[k, 2]),
                   d_dimer=float(at.get_distance(i, k, mic=True)), magmom=mag, removed_H=hd),
              open(f'runs/{tag}.json', 'w'), indent=1)
    print(tag, e, 'dz', at.positions[i, 2] - at.positions[k, 2])
