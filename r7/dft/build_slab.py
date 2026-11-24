#!/usr/bin/env python3
"""Build Si(100)-2x1:H slabs with one half-passivated dimer.

build(nx, ny, layers, side) -> Atoms
  nx, ny : surface cell in 1x1 units (nx even; dimers lie along x)
  layers : number of Si layers
  side   : 'L' or 'R' = which atom of the target dimer keeps its H
Bottom Si layer is dihydride-terminated; bottom two Si layers and their H are fixed.
"""
import sys
import numpy as np
from ase import Atoms
from ase.build import diamond100
from ase.constraints import FixAtoms
from ase.io import write

A0 = 5.431          # Si lattice constant (Å); replace with the PBE value from calibrate.py
D_SISI = 2.40       # starting monohydride dimer length (Å)
D_SIH = 1.50


def _neighbors(slab, i, rc=2.6):
    d = slab.get_distances(i, range(len(slab)), mic=True, vector=True)
    r = np.linalg.norm(d, axis=1)
    return [(j, d[j]) for j in np.where((r > 0.1) & (r < rc))[0]]


def build(nx=4, ny=4, layers=6, side='L', vacuum=15.0, a0=A0):
    s = diamond100('Si', (nx, ny, layers + 1), a=a0, vacuum=vacuum, periodic=True)
    s.pbc = (True, True, False)
    tags = s.get_tags()                       # 1 = top layer ... layers+1 = sacrificial bottom
    top = np.where(tags == 1)[0]
    # dimer axis: top atoms' bonds go to layer 2 along one in-plane axis; dimers form along the other
    j, v = [(j, v) for j, v in _neighbors(s, top[0]) if tags[j] == 2][0]
    bond_axis = int(np.argmax(np.abs(v[:2])))
    dim_axis = 1 - bond_axis
    pos = s.get_positions()
    # pair top atoms along dim_axis into dimers
    cell = s.cell.lengths()
    step = cell[dim_axis]/(nx if dim_axis == 0 else ny)
    dimers = []
    for i in top:
        k = int(round((pos[i, dim_axis] % cell[dim_axis])/step)) % (nx if dim_axis == 0 else ny)
        if k % 2 == 0:
            # partner: same other coordinate, next site along dim_axis
            target = pos[i].copy(); target[dim_axis] += step
            dist = [np.linalg.norm(((pos[m] - target + 0.5*np.append(cell[:2], 1e9)) % np.append(cell[:2], 1e9)) - 0.5*np.append(cell[:2], 1e9)) for m in top]
            dimers.append((i, top[int(np.argmin(dist))]))
    shift = (step - D_SISI)/2
    hpos = []
    for i, k in dimers:
        u = np.zeros(3); u[dim_axis] = 1.0
        pos[i] += shift*u; pos[k] -= shift*u
        for m, sgn in ((i, -1), (k, +1)):
            dvec = np.array([0, 0, 1.0]) + 0.38*sgn*u
            hpos.append((m, pos[m] + D_SIH*dvec/np.linalg.norm(dvec)))
    s.set_positions(pos)
    # bottom dihydride: replace bonds to sacrificial layer by H
    bottom = np.where(tags == layers)[0]
    sac = set(np.where(tags == layers + 1)[0])
    for i in bottom:
        for j, v in _neighbors(s, i):
            if j in sac:
                hpos.append((i, s.positions[i] + 1.48*v/np.linalg.norm(v)))
    sac_sorted = sorted(sac)
    remap = lambda old: old - int(np.searchsorted(sac_sorted, old))
    del s[sac_sorted]
    dimers = [(remap(i), remap(k)) for i, k in dimers]
    hpos = [(remap(o), p) for o, p in hpos]
    pos = s.get_positions()
    # target dimer: the one nearest the cell centre; remove one H to make a single DB
    c = s.cell.lengths()[:2]/2
    cen = [np.linalg.norm(((pos[i, :2] + pos[k, :2])/2) - c) for i, k in dimers]
    ti, tk = dimers[int(np.argmin(cen))]
    keep_on = ti if side == 'L' else tk
    db_on = tk if side == 'L' else ti
    for owner, p in hpos:
        if owner == db_on:
            continue
        s.append('H'); s.positions[-1] = p
    s.center(vacuum=vacuum, axis=2)
    tags = list(s.get_tags())
    zSi = s.positions[[a.index for a in s if a.symbol == 'Si'], 2]
    zcut = np.sort(np.unique(np.round(zSi, 1)))[1] + 0.3
    fixed = [a.index for a in s if a.position[2] < zcut]
    s.set_constraint(FixAtoms(indices=fixed))
    s.info.update(dict(target_dimer=(int(ti), int(tk)), H_on=int(keep_on), DB_on=int(db_on),
                       dim_axis=dim_axis, layers=layers, nx=nx, ny=ny))
    return s


if __name__ == '__main__':
    nx, ny, L = (int(x) for x in sys.argv[1:4]) if len(sys.argv) > 3 else (4, 4, 6)
    for side in 'LR':
        at = build(nx, ny, L, side)
        write(f'slab_{nx}x{ny}x{L}_{side}.extxyz', at)
        print(side, len(at), at.get_chemical_formula(), at.info)
