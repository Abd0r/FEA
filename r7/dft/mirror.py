"""Map an L structure to R by the mirror plane through the target dimer centre (perpendicular to the dimer axis).
Usage: python mirror.py IN.extxyz OUT.extxyz [REF.extxyz]  -> reorders atoms to match REF and prints max mismatch."""
import sys, numpy as np
from ase.io import read, write
from scipy.optimize import linear_sum_assignment
a = read(sys.argv[1]); ref = read(sys.argv[3]) if len(sys.argv) > 3 else None
ax = a.info['dim_axis']; i, k = a.info['target_dimer']
src = ref if ref is not None else a   # plane from the unrelaxed symmetric slab when given
yc = 0.5*(src.positions[i, ax] + src.positions[k, ax])
p = a.positions.copy(); p[:, ax] = 2*yc - p[:, ax]
b = a.copy(); b.positions = p; b.wrap()
if ref is not None:
    L = ref.cell.lengths()
    def d(x, y):
        v = x[:, None, :] - y[None, :, :]; v[..., :2] -= np.round(v[..., :2]/L[:2])*L[:2]; return np.linalg.norm(v, axis=2)
    out = b.copy()
    for sym in set(ref.get_chemical_symbols()):
        ir = [n for n, s in enumerate(ref.get_chemical_symbols()) if s == sym]
        ib = [n for n, s in enumerate(b.get_chemical_symbols()) if s == sym]
        D = d(ref.positions[ir], b.positions[ib]); r, c = linear_sum_assignment(D)
        out.positions[ir] = b.positions[np.array(ib)[c]]
        print(sym, 'max mismatch %.3f A' % D[r, c].max())
    out.set_constraint(ref.constraints); out.info = dict(ref.info); b = out
write(sys.argv[2], b)
