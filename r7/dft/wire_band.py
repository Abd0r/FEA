"""Band width of a bare-dimer wire (the FEA wire) on H-Si(100): one full dimer row depassivated.
The 4x4 cell holds 4 dimers along the row (x), so bands are folded 4x; the wire band width is the
energy span of the gap states localised on the row. 1D tight binding: W = 4 t.
Start: alternately buckled (+-0.2 A) so the relaxation can find a buckled ground state if one exists.
Usage: gpaw python wire_band.py   -> runs/wire_band.json"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW
from gpaw.mpi import world
from relax import make_calc
at = read('slab_4x4x6_L.extxyz')
i, k = at.info['target_dimer']; ax = int(at.info['dim_axis'])
z = at.positions[:, 2]; L = at.cell[ax, ax]
c = 0.5*(at.positions[i, ax] + at.positions[k, ax])
top = [n for n in range(len(at)) if at[n].symbol == 'Si' and z[n] > z[i] - 0.5]
row = [n for n in top if abs(((at.positions[n, ax] - c + L/2) % L) - L/2) < 2.0]
hs = [n for n, s in enumerate(at.get_chemical_symbols()) if s == 'H']
hrow = [n for n in hs if min(at.get_distance(m, n, mic=True) for m in row) < 1.8]
assert len(row) == 8 and len(hrow) == 7, (row, hrow)   # the L slab already lacks the H on atom 93
del at[hrow]                                  # H follow Si, so row indices are unchanged
# alternate buckling start: within each dimer one atom up, one down; alternate along the row
pairs = sorted([(m, n) for m in row for n in row if m < n and abs(at.positions[m, 0] - at.positions[n, 0]) < 0.5],
               key=lambda p: at.positions[p[0], 0])
for j, (m, n) in enumerate(pairs):
    s = 0.2*(-1)**j
    at.positions[m, 2] += s; at.positions[n, 2] -= s
at.calc = make_calc(at, 'dzp', 0, txt='runs/wire_band_relax.txt', kpts=(4, 1, 1))
opt = BFGS(at, trajectory='runs/wire_band.traj', logfile='runs/wire_band.log')
opt.run(fmax=0.03, steps=300)
e0 = at.get_potential_energy(); ef = at.calc.get_fermi_level()
write('runs/wire_band.extxyz', at)
buck = [float(at.positions[m, 2] - at.positions[n, 2]) for m, n in pairs]
# non-self-consistent bands along the row
ks = [(x, 0, 0) for x in np.linspace(0, 0.5, 9)]
bc = at.calc.fixed_density(kpts=ks, symmetry='off', txt='runs/wire_band_bands.txt')
nb = bc.get_number_of_bands()
# weight of each state on the row region (|psi|^2 within 3 A of the row centre in y, above the top layer - 1 A)
gd_shape = bc.get_pseudo_wave_function(band=0, kpt=0).shape
Y = np.linspace(0, at.cell[1, 1], gd_shape[1], endpoint=False)
Z = np.linspace(0, at.cell[2, 2], gd_shape[2], endpoint=False)
yc = at.positions[row, 1].mean(); ztop = at.positions[row, 2].min() - 1.0
dy = np.abs(((Y - yc + at.cell[1, 1]/2) % at.cell[1, 1]) - at.cell[1, 1]/2)
mask = (dy[None, :, None] < 3.0) & (Z[None, None, :] > ztop)
out = dict(E0=float(e0), Ef=float(ef), buckling=buck, kx=[k_[0] for k_ in ks], bands=[])
for kk in range(len(ks)):
    eps = bc.get_eigenvalues(kpt=kk, spin=0)
    for n in range(nb):
        if abs(eps[n] - ef) < 1.5:
            psi = bc.get_pseudo_wave_function(band=n, kpt=kk, spin=0)
            w = float((np.abs(psi)**2*mask).sum()/(np.abs(psi)**2).sum())
            out['bands'].append(dict(k=kk, n=n, e=float(eps[n] - ef), w_row=w))
occ = [b_['e'] for b_ in out['bands'] if b_['w_row'] > 0.4 and b_['e'] < 0]
emp = [b_['e'] for b_ in out['bands'] if b_['w_row'] > 0.4 and b_['e'] > 0]
if occ: out['pi_width_eV'] = max(occ) - min(occ); out['t_pi_eV'] = out['pi_width_eV']/4
if emp: out['pistar_width_eV'] = max(emp) - min(emp); out['t_pistar_eV'] = out['pistar_width_eV']/4
if world.rank == 0:
    json.dump(out, open('runs/wire_band.json', 'w'), indent=1)
    print('buckling', buck, 'pi*', out.get('pistar_width_eV'), 'pi', out.get('pi_width_eV'))
