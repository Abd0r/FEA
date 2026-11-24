import time, sys
from ase.io import read
from gpaw import GPAW
from gpaw.mpi import world
basis = sys.argv[1]
at = read('slab_4x4x6_L.extxyz')
calc = GPAW(mode='lcao', basis=basis, xc='PBE', kpts=(1, 1, 1), h=0.2,
            poissonsolver={'dipolelayer': 'xy'}, spinpol=True,
            occupations={'name': 'fermi-dirac', 'width': 0.05},
            txt=f'timing_{basis}.txt')
at.calc = calc
t = time.time(); e = at.get_potential_energy(); t1 = time.time() - t
f = at.get_forces()
if world.rank == 0:
    print(basis, 'E', round(e, 4), 'SCF+forces s', round(time.time() - t, 1), 'energy s', round(t1, 1), flush=True)
