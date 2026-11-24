"""CI-NEB for the intra-dimer H hop. Usage: gpaw python neb.py L.extxyz R.extxyz OUT [basis] [charge] [nimg]
Each image keeps its own calculator (so SCF restarts from its previous density). IDPP start, climb from the start."""
import sys, json, time
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE
from gpaw.mpi import world
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from relax import make_calc
iL, iR, out = sys.argv[1:4]
basis = sys.argv[4] if len(sys.argv) > 4 else 'dzp'
charge = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
n = int(sys.argv[6]) if len(sys.argv) > 6 else 3
a, b = read(iL), read(iR)
imgs = [a] + [b.copy() for _ in range(n)] + [b]
neb = NEB(imgs, climb=True, k=0.5)
init = sys.argv[7] if len(sys.argv) > 7 else None
if init:   # start interior images from an earlier converged band (e.g. the neutral one)
    band = read(init, index=':')
    for im, src in zip(imgs[1:-1], band[1:-1]):
        im.positions = src.positions
else:
    neb.interpolate('idpp', mic=True, apply_constraint=True)
from ase.calculators.singlepoint import SinglePointCalculator
eL = a.calc.results['energy'] if a.calc else None   # relaxed endpoints carry their energy from extxyz
for j, im in enumerate(imgs):
    if j in (0, n + 1) and eL is not None:
        e = im.calc.results['energy'] if im.calc else eL   # mirrored R: equal to L by exact symmetry
        im.calc = SinglePointCalculator(im, energy=e, forces=0*im.positions)
    else:
        im.calc = make_calc(im, basis, charge, txt=f'{out}_img{j}.txt')
t = time.time()
opt = FIRE(neb, trajectory=out + '.traj', logfile=out + '.log', maxstep=0.1)
opt.run(fmax=0.05, steps=int(__import__("os").environ.get("NEB_STEPS", 150)))
E = [im.get_potential_energy() for im in imgs]
conv = bool(opt.converged(neb.get_forces()))   # on every rank, before the rank-0 block
write(out + '_band.extxyz', imgs)
if world.rank == 0:
    json.dump(dict(energies=E, Ea_fwd=max(E) - E[0], Ea_bwd=max(E) - E[-1], dE=E[-1] - E[0],
                   steps=opt.nsteps, converged=conv, wall_s=time.time() - t, basis=basis, charge=charge),
              open(out + '.json', 'w'), indent=1)
