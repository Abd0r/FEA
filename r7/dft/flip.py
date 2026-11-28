"""Buckling flip of a bare dimer, and where its charge sits.

A bare dimer buckles (DFT: 0.73 A neutral, 0.42 A with one extra electron), so its charge sits nearer one atom.
Two numbers decide whether the side-tolerant Block survives this:
 1. offset: how far the charge centroid sits from the dimer centre, in units of the 0.1125 nm half-flip used
    by opt/side_flip.py (that script needs off <= 0.3 for the absorber budget). Hirshfeld charges on the relaxed
    buckled geometry; the stored electron is the difference rho(-1) - rho(0).
 2. barrier: CI-NEB between the buckled state and its mirror image (exact by the p(2x1) mirror plane through the
    dimer centre, so dE = 0 and the climbing image is the saddle).
The jellium is pinned to one z window for every image, as in htst.py, so energies are comparable.

Usage: gpaw python flip.py charges            -> runs/flip_charges.json
       gpaw python flip.py neb CHARGE(0|-1) [nimg]  -> runs/flip_neb_{q0|qm}.json
"""
import sys, os, json, time, gc
import numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE
from gpaw.mpi import world
from relax import make_calc

HERE = os.path.dirname(os.path.abspath(__file__))
REF = read(os.path.join(HERE, 'slab_4x4x6_L.extxyz'))
I, K = (int(x) for x in REF.info['target_dimer'])
AX = int(REF.info['dim_axis'])
YC = 0.5*(REF.positions[I, AX] + REF.positions[K, AX])      # mirror plane: unrelaxed dimer centre
HALF = 1.125                                                # A, half the flip distance in side_flip.py
MIX = {'backend': 'pulay', 'method': 'sum', 'beta': 0.02, 'nmaxold': 8, 'weight': 50}   # msr1 default oscillates near the symmetric image


def state(q):
    return read(os.path.join(HERE, 'runs', f"bare_buck_{'q0' if q == 0 else 'qm'}.extxyz"))


def zjell(at):
    return (at.positions[:, 2].max() + 5.0, at.cell[2, 2] - 2.0)


def magmoms(at, q):
    m = np.zeros(len(at))
    if q != 0:
        m[I] = m[K] = 0.5
    else:
        m[I], m[K] = 0.5, -0.5     # broken-symmetry seed: lets the symmetric (saddle) dimer go open-shell if lower;
                                   # a buckled closed-shell image collapses to M = 0 (zero seeds trapped the H hop, 27 Nov)
    return m


def mirror(a):
    """mirror image of a through the plane AX = YC, atoms reordered to match a; returns (b, max mismatch of fixed atoms)"""
    b = a.copy(); p = b.positions.copy(); p[:, AX] = 2*YC - p[:, AX]; b.positions = p; b.wrap()
    L = a.cell.lengths()
    out = a.copy()
    for sym in set(a.get_chemical_symbols()):
        ia = [n for n, s in enumerate(a.get_chemical_symbols()) if s == sym]
        v = a.positions[ia][:, None, :] - b.positions[ia][None, :, :]
        v[..., :2] -= np.round(v[..., :2]/L[:2])*L[:2]
        D = np.linalg.norm(v, axis=2); r, c = linear_sum_assignment(D)
        out.positions[np.array(ia)[r]] = b.positions[np.array(ia)[c]]
    fixed = a.constraints[0].index
    v = out.positions[fixed] - a.positions[fixed]; v[:, :2] -= np.round(v[:, :2]/L[:2])*L[:2]
    # the mirrored dimer must be the swap of the original one
    sw = abs(out.positions[I, 2] - a.positions[K, 2]) + abs(out.positions[K, 2] - a.positions[I, 2])
    return out, float(np.linalg.norm(v, axis=1).max()), float(sw)


def hirshfeld(at, q, zj, tag):
    from gpaw.analyse.hirshfeld import HirshfeldPartitioning
    at = at.copy(); at.set_initial_magnetic_moments(magmoms(at, q))
    at.calc = make_calc(at, 'dzp', q, txt=f'runs/flip_{tag}.txt', zjell=zj, legacy_gpaw=True)
    e = at.get_potential_energy()
    ch = HirshfeldPartitioning(at.calc).get_charges()
    return e, np.array(ch)


def near_mask(geo, r=4.0, depth=4.0):
    c = 0.5*(geo.positions[I, :2] + geo.positions[K, :2])
    m = np.linalg.norm(geo.positions[:, :2] - c, axis=1) < r
    return m & (geo.positions[:, 2] > geo.positions[:, 2].max() - depth)


def charges():
    """Site dipole along the dimer axis, p = sum_near q_a (y_a - y_c), for each charge state at its own relaxed
    geometry. A flip reverses p, i.e. changes it by 2p; side_flip.py moves a charge -e by 2*off*1.125 A, so the
    equivalent offset is |p|/(e*1.125 A). Also the vertical (fixed-geometry) electron centroid rho(-1)-rho(0)."""
    zj = zjell(state(-1))                                     # one jellium window for every calculation
    runs = {}
    for name, geo, q in (('qm_at_qm', state(-1), -1), ('q0_at_qm', state(-1), 0), ('q0_at_q0', state(0), 0)):
        e, ch = hirshfeld(geo, q, zj, 'hirsh_' + name)
        y = geo.positions[:, AX] - YC; m = near_mask(geo)
        runs[name] = dict(E=float(e), q_I=float(ch[I]), q_K=float(ch[K]), y_I=float(y[I]), y_K=float(y[K]),
                          dz=float(geo.positions[I, 2] - geo.positions[K, 2]),
                          p_near_eA=float(np.sum(ch[m]*y[m])), p_dimer_eA=float(ch[I]*y[I] + ch[K]*y[K]),
                          q_near=float(ch[m].sum()), ch=ch, y=y, m=m)
    a, b = runs['qm_at_qm'], runs['q0_at_qm']
    dq = a['ch'] - b['ch']; m = a['m']
    cen = float(np.sum(dq[m]*a['y'][m])/np.sum(dq[m]))
    out = dict(jellium=list(zj), half_flip_A=HALF,
               offset_occupied=abs(a['p_near_eA'])/HALF, offset_empty=abs(runs['q0_at_q0']['p_near_eA'])/HALF,
               electron_centroid_A=cen, electron_offset_vertical=abs(cen)/HALF,
               electron_in_near=float(dq[m].sum()), electron_total=float(dq.sum()),
               note='Hirshfeld; atomic dipoles neglected. offset = |site dipole| / (e * 1.125 A); side_flip needs <= 0.3',
               states={k: {kk: vv for kk, vv in v.items() if kk not in ('ch', 'y', 'm')} for k, v in runs.items()})
    if world.rank == 0:
        json.dump(out, open('runs/flip_charges.json', 'w'), indent=1)
        print(json.dumps(out, indent=1))


def neb(q, n):
    tag = f"flip_neb_{'q0' if q == 0 else 'qm'}"
    a = state(q)
    b, mis, sw = mirror(a)
    if world.rank == 0:
        print(f'mirror: fixed-atom mismatch {mis:.2e} A, dimer swap mismatch {sw:.2e} A', flush=True)
    assert mis < 1e-3 and sw < 1e-3, (mis, sw)
    zj = zjell(a)
    imgs = [a] + [a.copy() for _ in range(n)] + [b]
    band = NEB(imgs, climb=True, k=0.5)
    band.interpolate('idpp', mic=True, apply_constraint=True)
    tr = f'runs/{tag}.traj'
    if os.path.exists(tr) and os.path.getsize(tr) > 0:              # restart from the last saved band
        last = read(tr, index=f'-{n + 2}:')
        if len(last) == n + 2:
            for im, old in zip(imgs[1:-1], last[1:-1]):
                im.positions = old.positions
            world.barrier()
            if world.rank == 0:
                os.rename(tr, tr.replace('.traj', f'_part{int(time.time())}.traj')); print(tag, 'restarted', flush=True)
            world.barrier()
    # endpoints one at a time on all ranks, then freed; cached for restarts (28 Nov: five resident calculators
    # needed ~15-22 GB on the 16 GB machine)
    ends = f'runs/{tag}_ends.json'
    if os.path.exists(ends):
        E_end = json.load(open(ends))
    else:
        E_end = []
        for j in (0, n + 1):
            x = imgs[j].copy(); x.set_initial_magnetic_moments(magmoms(x, q))
            x.calc = make_calc(x, 'dzp', q, txt=f'runs/{tag}_img{j}.txt', zjell=zj, mixer=MIX, maxiter=500)
            E_end.append(float(x.get_potential_energy())); x.calc = None; del x; gc.collect()
        world.barrier()
        if world.rank == 0:
            json.dump(E_end, open(ends, 'w'))
    if world.rank == 0:
        print('endpoint energies', E_end, 'difference (symmetry check)', E_end[1] - E_end[0], flush=True)
    from ase.calculators.singlepoint import SinglePointCalculator
    for j in (0, n + 1):
        imgs[j].calc = SinglePointCalculator(imgs[j], energy=E_end[0 if j == 0 else 1], forces=0*imgs[j].positions)
    # moving images in parallel: world.size/n ranks per image
    assert world.size % n == 0, f'use a multiple of {n} MPI ranks'
    per = world.size // n; mine = world.rank // per
    comm = world.new_communicator(list(range(mine*per, (mine + 1)*per)))
    for j in range(1, n + 1):
        imgs[j].set_initial_magnetic_moments(magmoms(imgs[j], q))
        if j - 1 == mine:
            imgs[j].calc = make_calc(imgs[j], 'dzp', q, txt=f'runs/{tag}_img{j}.txt', zjell=zj, mixer=MIX,
                                     maxiter=500, communicator=comm)
    band = NEB(imgs, climb=True, k=0.5, parallel=True)
    t = time.time()
    opt = FIRE(band, trajectory=f'runs/{tag}.traj', logfile=f'runs/{tag}.log', maxstep=0.1)
    opt.run(fmax=0.05, steps=int(os.environ.get('NEB_STEPS', 200)))
    conv = bool(opt.converged(band.get_forces()))
    E = [float(e) for e in band.energies]
    dz = [float(im.positions[I, 2] - im.positions[K, 2]) for im in imgs]
    write(f'runs/{tag}_band.extxyz', imgs)
    if world.rank == 0:
        json.dump(dict(charge=q, n_images=n, energies=E, Ea=max(E) - E[0], Ea_from_mean_end=max(E) - 0.5*sum(E_end),
                       dE_ends=E_end[1] - E_end[0], dz_path=dz, steps=opt.nsteps, converged=conv,
                       mirror_mismatch_A=mis, wall_s=time.time() - t, jellium=list(zj)),
                  open(f'runs/{tag}.json', 'w'), indent=1)
        print(tag, 'Ea', max(E) - E[0], 'converged', conv)


if __name__ == '__main__':
    if sys.argv[1] == 'charges':
        charges()
    elif sys.argv[1] == 'mirrorcheck':
        for q in (0, -1):
            b, mis, sw = mirror(state(q)); print(q, 'fixed mismatch', mis, 'swap', sw)
    else:
        neb(float(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 3)
