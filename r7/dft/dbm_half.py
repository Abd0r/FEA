"""DB- intra-dimer hop: is the 1.43 eV saddle the lowest?

The 5-image DB- band (neb_qm5) has mirror-image endpoints, but its two halves follow different paths
(images j and mirror(6-j) differ by up to 0.33 A): the left half climbed to 1.43 eV, the right half was not
climbing and its highest image sits at 1.33 eV. By the mirror symmetry the right-half path, mirrored, is also a
path from L to the bridged intermediate. This script
 1. relaxes the bridged intermediate M (image 3) as a minimum;
 2. runs a climbing-image NEB from L to M started from the mirrored right half;
 3. recomputes L and the old saddle (image 2) with the same pinned jellium window,
so the DB- barrier is min(old saddle, new saddle), all on one energy reference.
Usage: gpaw python dbm_half.py   -> runs/dbm_half.json
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import linear_sum_assignment
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
YC = 0.5*(REF.positions[I, AX] + REF.positions[K, AX])
MIX = {'backend': 'pulay', 'method': 'sum', 'beta': 0.02, 'nmaxold': 8, 'weight': 50}
Q = -1


def mirror_pos(p):
    p = p.copy(); p[:, AX] = 2*YC - p[:, AX]; return p


def sigma(L, R):
    """atom permutation of the mirror: R[sigma[a]] ~ mirror(L[a])"""
    m = mirror_pos(L.positions); cell = L.cell.lengths(); s = np.arange(len(L))
    for sym in set(L.get_chemical_symbols()):
        idx = np.array([n for n, x in enumerate(L.get_chemical_symbols()) if x == sym])
        v = m[idx][:, None] - R.positions[idx][None]; v[..., :2] -= np.round(v[..., :2]/cell[:2])*cell[:2]
        D = np.linalg.norm(v, axis=2); r, c = linear_sum_assignment(D)
        assert D[r, c].max() < 0.05, ('mirror map', sym, D[r, c].max())
        s[idx[r]] = idx[c]
    return s


def mirrored(X, s):
    Y = X.copy(); p = mirror_pos(X.positions); Y.positions[s] = p; Y.wrap(); return Y


def calc(at, txt, zj):
    m = np.zeros(len(at)); m[I] = m[K] = 0.5; at.set_initial_magnetic_moments(m)
    return make_calc(at, 'dzp', Q, txt=txt, zjell=zj, mixer=MIX, maxiter=500)


band = read(os.path.join(HERE, 'runs', 'neb_qm5_band.extxyz'), index=':')
L, R = band[0], band[-1]
s = sigma(L, R)
zj = (L.positions[:, 2].max() + 5.0, L.cell[2, 2] - 2.0)
out = dict(jellium=list(zj))
# asymmetry of the old band
cell = L.cell.lengths()
def dev(A, B):
    v = A.positions - B.positions; v[:, :2] -= np.round(v[:, :2]/cell[:2])*cell[:2]; return float(np.linalg.norm(v, axis=1).max())
out['old_band_asymmetry_A'] = [dev(band[j], mirrored(band[6 - j], s)) for j in range(7)]

# 1. bridged intermediate
M = band[3].copy(); M.calc = calc(M, 'runs/dbm_bridge.txt', zj)
opt = BFGS(M, trajectory='runs/dbm_bridge.traj', logfile='runs/dbm_bridge.log'); opt.run(fmax=0.03, steps=300)
E_M = M.get_potential_energy(); write('runs/dbm_bridge.extxyz', M)
out.update(E_bridge=float(E_M), bridge_steps=opt.nsteps, bridge_asymmetry_A=dev(M, mirrored(M, s)))

# 3. same-reference energies of L and the old saddle
for name, X in (('L', L), ('old_saddle', band[2])):
    Y = X.copy(); Y.calc = calc(Y, f'runs/dbm_sp_{name}.txt', zj); out[f'E_{name}'] = float(Y.get_potential_energy())

# 2. CI-NEB L -> M from the mirrored right half
m5, m4 = mirrored(band[5], s), mirrored(band[4], s)
mid = m5.copy(); mid.positions = 0.5*(m5.positions + m4.positions)
imgs = [L.copy(), m5, mid, m4, M.copy()]
for j, im in enumerate(imgs):
    im.calc = calc(im, f'runs/dbm_half_img{j}.txt', zj)
from ase.calculators.singlepoint import SinglePointCalculator
imgs[0].calc = SinglePointCalculator(imgs[0], energy=out['E_L'], forces=0*imgs[0].positions)
imgs[-1].calc = SinglePointCalculator(imgs[-1], energy=E_M, forces=0*imgs[-1].positions)
nb = NEB(imgs, climb=True, k=0.5, method='improvedtangent')
t0 = time.time()
opt = FIRE(nb, trajectory='runs/dbm_half.traj', logfile='runs/dbm_half.log', maxstep=0.1)
opt.run(fmax=0.05, steps=int(os.environ.get('NEB_STEPS', 200)))
E = [im.get_potential_energy() for im in imgs]
write('runs/dbm_half_band.extxyz', imgs)
out.update(new_path_energies=[e - out['E_L'] for e in E], new_saddle=max(E) - out['E_L'],
           old_saddle=out['E_old_saddle'] - out['E_L'], bridge=E_M - out['E_L'],
           neb_steps=opt.nsteps, converged=bool(opt.converged(nb.get_forces())), wall_s=time.time() - t0)
out['Ea_DBminus'] = min(out['new_saddle'], out['old_saddle'])
if world.rank == 0:
    json.dump(out, open('runs/dbm_half.json', 'w'), indent=1)
    print('old saddle %.3f  new saddle %.3f  bridge %.3f  converged %s' % (out['old_saddle'], out['new_saddle'],
                                                                         out['bridge'], out['converged']))
