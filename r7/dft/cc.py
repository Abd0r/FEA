"""Configuration-coordinate single points for optimisation 1 (lattice-relaxation barrier).
E(q @ geometry of q') for q, q' in {0,-1,+1}, all charged runs using ONE fixed jellium slab (from q0_L)
so energies at the same charge but different geometry are directly comparable.
Relaxation (reorganisation) energies: lam(q at q' geom) = E(q @ q') - E(q @ q)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read
from gpaw.mpi import world
from relax import make_calc
geo = {0: read('runs/q0_L.extxyz'), -1: read('runs/qm_L.extxyz'), 1: read('runs/qp_L.extxyz')}
ref = geo[0]; zj = (ref.positions[:, 2].max() + 5.0, ref.cell[2, 2] - 2.0)
jobs = [(0, -1), (0, 1), (-1, 0), (-1, -1), (1, 0), (1, 1)]   # (charge, geometry)
E = {}
for q, g in jobs:
    at = geo[g].copy(); at.calc = make_calc(at, 'dzp', q, txt=f'runs/cc_q{q}_g{g}.txt', zjell=zj)
    E[f'{q}@{g}'] = at.get_potential_energy()
E['0@0'] = read('runs/q0_L.extxyz').get_potential_energy()
lam = {'neutral, at DB- geometry': E['0@-1'] - E['0@0'], 'DB- at neutral geometry': E['-1@0'] - E['-1@-1'],
       'neutral, at DB+ geometry': E['0@1'] - E['0@0'], 'DB+ at neutral geometry': E['1@0'] - E['1@1']}
if world.rank == 0:
    json.dump(dict(E=E, lam=lam, zjell=zj), open('runs/cc.json', 'w'), indent=1)
