"""Convergence of the intra-dimer barrier against the three approximations a computational referee will question.

Single points on the existing geometries (endpoint and climbing image of the CI-NEB bands):
  neutral:  Gamma LCAO-dzp (reference), 2x2x1 k-points (LCAO-dzp), plane waves 400 eV at Gamma;
  charged:  jellium window as used, and moved 3 A further from the slab (field / image artefacts).
Each barrier is E(saddle) - E(endpoint) evaluated with identical settings.
Usage: gpaw python conv.py [case ...]   -> runs/conv.json (cases: lcao kpt pw jell)
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ase.io import read
from gpaw.mpi import world
from relax import make_calc

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda f, i=None: read(os.path.join(HERE, 'runs', f), index=i) if i is not None else read(os.path.join(HERE, 'runs', f))
MIX = {'backend': 'pulay', 'method': 'sum', 'beta': 0.02, 'nmaxold': 8, 'weight': 50}
OUT = os.path.join(HERE, 'runs', 'conv.json')


def energy(at, q, tag, **kw):
    at = at.copy()
    ref = read(os.path.join(HERE, 'slab_4x4x6_L.extxyz')); i, k = (int(x) for x in ref.info['target_dimer'])
    m = np.zeros(len(at)); m[i] = m[k] = 0.5; at.set_initial_magnetic_moments(m)
    at.calc = make_calc(at, kw.pop('basis', 'dzp'), q, txt=f'runs/conv_{tag}.txt', mixer=MIX, maxiter=500, **kw)
    t = time.time(); e = at.get_potential_energy()
    return float(e), time.time() - t


def barrier(endpoint, saddle, q, tag, **kw):
    e0, t0 = energy(endpoint, q, tag + '_end', **dict(kw)); e1, t1 = energy(saddle, q, tag + '_sad', **dict(kw))
    return dict(Ea=e1 - e0, E_end=e0, E_saddle=e1, wall_s=t0 + t1)


q0_end, q0_sad = R('q0_L.extxyz'), R('neb_q0_band.extxyz', 2)
qm_end, qm_sad = R('qm_L.extxyz'), R('neb_qm5_band.extxyz', 2)
zj = (qm_end.positions[:, 2].max() + 5.0, qm_end.cell[2, 2] - 2.0)
cases = sys.argv[1:] or ['lcao', 'kpt', 'jell', 'pw']
out = json.load(open(OUT)) if os.path.exists(OUT) else {}
for c in cases:
    if c == 'lcao':
        out['q0_gamma_lcao'] = barrier(q0_end, q0_sad, 0, 'q0_lcao')
    elif c == 'kpt':
        out['q0_k221_lcao'] = barrier(q0_end, q0_sad, 0, 'q0_k221', kpts=(2, 2, 1))
    elif c == 'pw':
        out['q0_gamma_pw400'] = barrier(q0_end, q0_sad, 0, 'q0_pw', mode={'name': 'pw', 'ecut': 400})
    elif c == 'jell':
        out['qm_jell_as_used'] = barrier(qm_end, qm_sad, -1, 'qm_j0', zjell=zj)
        out['qm_jell_plus3A'] = barrier(qm_end, qm_sad, -1, 'qm_j3', zjell=(zj[0] + 3.0, zj[1]))
    if world.rank == 0:
        json.dump(out, open(OUT, 'w'), indent=1)
        print(c, {k: round(v['Ea'], 4) for k, v in out.items()}, flush=True)
