"""Redo the two NEUTRAL configuration-coordinate single points (neutral at the DB- and DB+ geometries)
with a starting moment of 1 muB on the DB atom. In cc.py they started from zero moment and converged to a
non-magnetic solution, the wrong state for a neutral DB. Merges the result into runs/cc.json as cc_fixed."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read
from gpaw.mpi import world
from relax import make_calc
cc = json.load(open('runs/cc.json'))
E = dict(cc['E'])
for g, f in ((-1, 'qm_L'), (1, 'qp_L')):
    at = read(f'runs/{f}.extxyz'); m = [0.0]*len(at); m[93] = 1.0; at.set_initial_magnetic_moments(m)
    at.calc = make_calc(at, 'dzp', 0, txt=f'runs/cc_fix_q0_g{g}.txt')
    E[f'0@{g}'] = at.get_potential_energy()
    E[f'0@{g}_magmom'] = float(at.calc.get_magnetic_moment())
lam = {'neutral, at DB- geometry': E['0@-1'] - E['0@0'], 'DB- at neutral geometry': E['-1@0'] - E['-1@-1'],
       'neutral, at DB+ geometry': E['0@1'] - E['0@0'], 'DB+ at neutral geometry': E['1@0'] - E['1@1']}
if world.rank == 0:
    cc.update(E_fixed=E, lam=lam, lam_unfixed=cc['lam'],
              note='neutral single points redone spin-polarized (cc_fix.py); lam_unfixed kept for record')
    json.dump(cc, open('runs/cc.json', 'w'), indent=1)
