"""Spin-polarized single points for the neutral band (endpoint + saddle), started with a moment on the target dimer.
The q0 relax/NEB started from zero moment and stayed spin-unpolarized. Correction:
Ea(spin) ~ Ea(q0 run) + [dE_spin(saddle) - dE_spin(min)], dE_spin = E(spin-pol) - E(q0 run)."""
import sys, json
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from ase.io import read
from gpaw.mpi import world
from relax import make_calc
band = read('runs/neb_q0_band.extxyz', index=':')
E0 = json.load(open('runs/neb_q0.json'))['energies']
isad = max(range(len(E0)), key=lambda j: E0[j])
out = {}
for j in (0, isad):
    at = band[j].copy(); i, k = at.info['target_dimer']
    m = [0.0]*len(at); m[i] = m[k] = 0.5; at.set_initial_magnetic_moments(m)
    at.calc = make_calc(at, 'dzp', 0, txt=f'runs/spin_img{j}.txt')
    e = at.get_potential_energy(); mag = at.calc.get_magnetic_moment()
    out[j] = dict(E_spin=e, E_q0=E0[j], dE_spin=e - E0[j], magmom=mag)
corr = out[isad]['dE_spin'] - out[0]['dE_spin']
if world.rank == 0:
    json.dump(dict(points=out, saddle_image=isad, Ea_q0=E0[isad] - E0[0], Ea_spin_corrected=E0[isad] - E0[0] + corr),
              open('runs/spin_fix.json', 'w'), indent=1)
