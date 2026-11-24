"""Add a neighbour DB on Si 92: the next dimer along the row, on the same side as target atom k=93.
Distances: 3.82 A to k, 4.44 A to i (a point charge there makes the k-side DB costlier for like charges).
Writes tilt_{L,R}.extxyz from the relaxed DB- endpoints. L = H on i (DB on k, NEAR the neighbour); R = H on k (DB on i, FAR)."""
import numpy as np
from ase.io import read, write
NB = 92
for s in 'LR':
    a = read(f'runs/qm_{s}.extxyz')
    h = [n for n, x in enumerate(a.get_chemical_symbols()) if x == 'H']
    d = a.get_distances(NB, h, mic=True); j = h[int(np.argmin(d))]
    assert d.min() < 1.7, d.min()
    a.calc = None
    del a[j]                                   # H indices > all Si, so Si indices (target 89/93, NB 92) are unchanged
    i, k = a.info['target_dimer']; db = k if s == 'L' else i
    m = np.zeros(len(a)); m[db] = 1.0; m[NB] = 1.0  # used for q=0 (two DB0); harmless otherwise
    a.set_initial_magnetic_moments(m); a.info['neighbour_DB'] = NB
    write(f'tilt_{s}.extxyz', a)
    print(s, a.get_chemical_formula(), 'target DB on', db, 'dist to NB %.2f' % a.get_distance(db, NB, mic=True))
