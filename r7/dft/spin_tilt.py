"""Neutral tilt, other spin state at each relaxed geometry (single points).
The relaxations ended in different states: near (L) paired, moment 0; far (R) triplet, moment 2.
Here: L as a triplet (moment fixed at 2), R as broken-symmetry singlet (+1/-1 start). The lower state of
each configuration is the one to compare."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read
from gpaw import GPAW
from gpaw.mpi import world
out = {}
# L triplet: converged in runs/spin_tilt_L.txt (Extrapolated -643.740426 eV, moment 2.000); GPAW then crashed while
# printing work functions (fixmagmom + dipole correction assertion), so its energy is taken from that log.
out['L'] = dict(E=-643.740426, magmom=2.000007, state='triplet', source='runs/spin_tilt_L.txt')
for s, mom, fix in (('R', (1.0, -1.0), False),):
    at = read(f'runs/tilt_q0_{s}.extxyz'); i, k = at.info['target_dimer']; db = k if s == 'L' else i
    m = [0.0]*len(at); m[db], m[92] = mom; at.set_initial_magnetic_moments(m)
    # same settings as relax.make_calc (neutral), with the total moment fixed for the triplet
    occ = {'name': 'fermi-dirac', 'width': 0.05, 'fixmagmom': True} if fix else {'name': 'fermi-dirac', 'width': 0.05}
    at.calc = GPAW(mode='lcao', basis='dzp', xc='PBE', kpts=(1, 1, 1), h=0.2, spinpol=True, occupations=occ,
                   poissonsolver={'dipolelayer': 'xy'}, txt=f'runs/spin_tilt_{s}.txt')
    out[s] = dict(E=at.get_potential_energy(), magmom=float(at.calc.get_magnetic_moment()), state='triplet' if fix else 'broken-symmetry')
if world.rank == 0:
    json.dump(out, open('runs/spin_tilt.json', 'w'), indent=1)
