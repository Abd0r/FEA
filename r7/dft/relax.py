"""Relax an endpoint. Usage: gpaw python relax.py IN.extxyz OUT_PREFIX [basis] [charge]
Charge != 0 uses a compensating jellium slab in the upper vacuum, plus a dipole correction."""
import sys, json, time
from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW
from gpaw.mpi import world


def make_calc(at, basis='dzp', charge=0, txt='-', zjell=None, **extra):
    kw = dict(mode='lcao', basis=basis, xc='PBE', kpts=(1, 1, 1), h=0.2, spinpol=True,
              occupations={'name': 'fermi-dirac', 'width': 0.05},
              poissonsolver={'dipolelayer': 'xy'}, txt=txt)
    kw.update(extra)
    if charge:
        from gpaw.jellium import JelliumSlab
        ztop = at.positions[:, 2].max()
        z1, z2 = zjell if zjell else (ztop + 5.0, at.cell[2, 2] - 2.0)
        # JelliumSlab(n) adds n electrons to the slab AND an equal positive background, so the cell stays neutral.
        # Do not also pass charge=: that double counts (net charge left -> Poisson solver refuses).
        kw.update(background_charge=JelliumSlab(-charge, z1=z1, z2=z2))
    return GPAW(**kw)


if __name__ == '__main__':
    inp, out = sys.argv[1], sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'dzp'
    charge = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    at = read(inp)
    at.calc = make_calc(at, basis, charge, txt=out + '.txt')
    t = time.time()
    opt = BFGS(at, trajectory=out + '.traj', logfile=out + '.log')
    opt.run(fmax=0.03, steps=300)
    # every rank evaluates everything first: a calculator call made by rank 0 alone deadlocks MPI
    e = at.get_potential_energy()
    f = at.get_forces()
    fmax = float((f**2).sum(axis=1).max()**0.5)
    write(out + '.extxyz', at)
    if world.rank == 0:
        json.dump(dict(input=inp, basis=basis, charge=charge, energy_eV=e,
                       fmax=fmax, steps=opt.nsteps,
                       wall_s=time.time() - t, converged=fmax < 0.03),
                  open(out + '.json', 'w'), indent=1)
