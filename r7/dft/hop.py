"""S3: an H atom from a neighbouring dimer (same row) hops INTO a bare dimer.

This is the event that damages a bare-dimer Block site: afterwards the site's dimer holds one H and the DB pair is split
between the site and its neighbour. The paper so far uses Hill's clean-surface intrarow barrier (1.75 eV); this
computes it on H-Si(100) for the bare-dimer site, neutral and charged, to both atoms of the buckled dimer.
The reverse barrier (final -> initial) says whether the damage heals itself (S5).

Targets (indices from bare_buck_*.extxyz): the up atom I=89 receives H104 from Si88; the down atom K=93 receives H105
from Si92. The other two neighbours (Si90, Si94) are mirror images along the row.
Usage: gpaw python hop.py relax CHARGE(0|-1) TARGET(up|down)   -> runs/hop_{q}_{t}_final.extxyz
       gpaw python hop.py neb   CHARGE(0|-1) TARGET(up|down) [nimg]  -> runs/hop_{q}_{t}.json
"""
import sys, os, json, time, gc
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS, FIRE
from gpaw.mpi import world
from relax import make_calc

HERE = os.path.dirname(os.path.abspath(__file__))
REF = read(os.path.join(HERE, 'slab_4x4x6_L.extxyz'))
I, K = (int(x) for x in REF.info['target_dimer'])
AX = int(REF.info['dim_axis'])
MIX = {'backend': 'pulay', 'method': 'sum', 'beta': 0.02, 'nmaxold': 8, 'weight': 50}
HOP = {'up': (I, 88, 104), 'down': (K, 92, 105)}      # (receiving Si, donor Si, moving H)
H_ON_I = 106                                          # index of the H on Si89 in slab_4x4x6_L


def extra(q):
    """Fixed total moment (fixmagmom) crashes this GPAW with the dipole layer (two Fermi levels, 28 Nov): free moment."""
    return dict(mixer=MIX, maxiter=500)


def tags(q, t):
    return f"hop_{'q0' if q == 0 else 'qm'}_{t}"


def initial(q):
    return read(os.path.join(HERE, 'runs', f"bare_buck_{'q0' if q == 0 else 'qm'}.extxyz"))


def zjell(at):
    return (at.positions[:, 2].max() + 5.0, at.cell[2, 2] - 2.0)


def magmoms(at, q, t):
    m = np.zeros(len(at))
    rec, don, _ = HOP[t]
    other = K if rec == I else I                     # the two DBs left after the hop: donor Si and the other dimer atom
    if q != 0:
        m[don] = 0.5; m[other] = 0.5
    else:
        # neutral: two singly occupied DBs. Seed a broken-symmetry singlet (+-1), continuous with the closed-shell
        # bare dimer at the start; a zero seed stays on the closed-shell surface, 0.34 eV above (found 27 Nov).
        # The triplet (what the final-state relax fell into) is recorded separately from hop_*_final.json.
        m[don] = 1.0; m[other] = -1.0
    return m


def guess_final(a, t):
    """H moved onto the receiving atom, along the Si-H bond direction of the hydrogenated surface."""
    rec, don, h = HOP[t]
    assert a.get_distance(don, h, mic=True) < 1.7, 'moving H is not on the donor'
    b = a.copy()
    p = REF.positions[H_ON_I] - REF.positions[I]                 # Si-H bond vector on Si89 in the H-terminated slab
    if rec == K:
        p = p.copy(); p[AX] = -p[AX]                             # mirror for the other dimer atom
    b.positions[h] = b.positions[rec] + p
    return b


def relax_final(q, t):
    a = initial(q); zj = zjell(a)
    b = guess_final(a, t)
    tag = tags(q, t)
    nf = f"runs/{tags(0, t)}_final.extxyz"
    if q != 0 and os.path.exists(nf):                              # charged: start from the relaxed neutral final state
        b.positions = read(nf).positions
    tr = f'runs/{tag}_final.traj'
    if os.path.exists(tr) and os.path.getsize(tr) > 0:              # restart from the last relaxed geometry
        b.positions = read(tr, index=-1).positions
        world.barrier()
        if world.rank == 0:
            os.rename(tr, tr.replace('.traj', f'_part{int(time.time())}.traj')); print(tag, 'relax restarted', flush=True)
        world.barrier()
    b.set_initial_magnetic_moments(magmoms(b, q, t))
    b.calc = make_calc(b, 'dzp', q, txt=f'runs/{tag}_final.txt', zjell=zj, **extra(q))
    opt = BFGS(b, trajectory=f'runs/{tag}_final.traj', logfile=f'runs/{tag}_final.log')
    opt.run(fmax=0.03, steps=300)
    e = b.get_potential_energy(); f = b.get_forces()
    write(f'runs/{tag}_final.extxyz', b)
    if world.rank == 0:
        json.dump(dict(energy_eV=e, fmax=float(np.linalg.norm(f, axis=1).max()), steps=opt.nsteps, jellium=list(zj)),
                  open(f'runs/{tag}_final.json', 'w'), indent=1)


def neb(q, t, n):
    tag = tags(q, t)
    a = initial(q); zj = zjell(a)
    b = read(f'runs/{tag}_final.extxyz')
    imgs = [a] + [a.copy() for _ in range(n)] + [b]
    band = NEB(imgs, climb=True, k=0.5, method='improvedtangent')
    band.interpolate('idpp', mic=True, apply_constraint=True)
    tr = f'runs/{tag}.traj'
    if os.path.exists(tr) and os.path.getsize(tr) > 0:              # restart from the last saved band (crash / disk full)
        last = read(tr, index=f'-{n + 2}:')
        if len(last) == n + 2:
            for im, old in zip(imgs[1:-1], last[1:-1]):
                im.positions = old.positions
            world.barrier()
            if world.rank == 0:
                os.rename(tr, tr.replace('.traj', f'_part{int(time.time())}.traj'))
                print(tag, 'restarted from', tr, flush=True)
            world.barrier()
    def mm(j, im):
        return magmoms(im, q, t) if j else (np.where(np.isin(np.arange(len(im)), [I, K]), 0.5, 0.0) if q else np.zeros(len(im)))
    # endpoints: on all ranks, one at a time, then freed (cached so a restart skips them)
    ends = f'runs/{tag}_ends.json'
    if os.path.exists(ends):
        E_end = json.load(open(ends))
    else:
        E_end = []
        for j in (0, n + 1):
            x = imgs[j].copy(); x.set_initial_magnetic_moments(mm(j, x))
            x.calc = make_calc(x, 'dzp', q, txt=f'runs/{tag}_img{j}.txt', zjell=zj, **extra(q))
            E_end.append(float(x.get_potential_energy())); x.calc = None; del x; gc.collect()
        world.barrier()
        if world.rank == 0:
            json.dump(E_end, open(ends, 'w'))
    E_triplet_final = json.load(open(f'runs/{tag}_final.json'))['energy_eV'] if q == 0 else None
    from ase.calculators.singlepoint import SinglePointCalculator
    for j, e in ((0, E_end[0]), (n + 1, E_end[1])):
        imgs[j].calc = SinglePointCalculator(imgs[j], energy=e, forces=0*imgs[j].positions)
    # moving images in parallel: world.size/n ranks per image, each rank holds ONE image (memory / n)
    assert world.size % n == 0, f'use a multiple of {n} MPI ranks'
    per = world.size // n; mine = world.rank // per
    comm = world.new_communicator(list(range(mine*per, (mine + 1)*per)))
    for j in range(1, n + 1):
        imgs[j].set_initial_magnetic_moments(mm(j, imgs[j]))
        if j - 1 == mine:
            imgs[j].calc = make_calc(imgs[j], 'dzp', q, txt=f'runs/{tag}_img{j}.txt', zjell=zj,
                                     communicator=comm, **extra(q))
    band = NEB(imgs, climb=True, k=0.5, method='improvedtangent', parallel=True)
    t0 = time.time()
    opt = FIRE(band, trajectory=f'runs/{tag}.traj', logfile=f'runs/{tag}.log', maxstep=0.1)
    opt.run(fmax=0.05, steps=int(os.environ.get('NEB_STEPS', 200)))
    conv = bool(opt.converged(band.get_forces()))
    E = [float(e) for e in band.energies]
    write(f'runs/{tag}_band.extxyz', imgs)
    if world.rank == 0:
        json.dump(dict(charge=q, target=t, n_images=n, energies=E, Ea_fwd=max(E) - E[0], Ea_bwd=max(E) - E[-1],
                       dE=E[-1] - E[0], E_final_relax=E_triplet_final, steps=opt.nsteps, converged=conv, wall_s=time.time() - t0, jellium=list(zj),
                       hop=dict(zip(('receiver', 'donor', 'H'), HOP[t]))),
                  open(f'runs/{tag}.json', 'w'), indent=1)
        print(tag, 'Ea_fwd', max(E) - E[0], 'Ea_bwd', max(E) - E[-1], 'converged', conv)


if __name__ == '__main__':
    mode, q, t = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    if mode == 'check':
        a = initial(q); b = guess_final(a, t); rec, don, h = HOP[t]
        print(t, 'H-receiver %.3f' % b.get_distance(h, rec, mic=True), 'H-donor %.3f' % b.get_distance(h, don, mic=True),
              'closest other atom %.3f' % min(b.get_distance(h, k, mic=True) for k in range(len(b)) if k not in (h, rec)))
    elif mode == 'relax':
        relax_final(q, t)
    else:
        neb(q, t, int(sys.argv[4]) if len(sys.argv) > 4 else 3)
