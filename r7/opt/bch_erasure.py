#!/usr/bin/env python3
"""Erasure capability of the FEA word code under maximum-likelihood (Gaussian-elimination) decoding.

The code is the narrow-sense binary BCH code of length 127 and designed distance 17 (t = 8, generator degree 56,
k = 71), shortened by 7 to (120, 64). A pattern of s erased positions is decodable iff no nonzero codeword is
supported inside it, i.e. iff the generator matrix restricted to the n - s surviving positions has full rank 64.
Patterns of s <= 16 always decode (minimum distance >= 17); this script measures P_fail(s) by Monte Carlo for
s > 16 and compares it with the random-code value 2^-(n-k-s) used in the chip model.
Usage: python3 bch_erasure.py [trials_per_s]  -> bch_erasure.json
"""
import sys, json, os, random

M, PRIM = 7, 0b10001001                  # GF(2^7), x^7 + x^3 + 1
N0, T = 127, 8
SHORT = 7

def gf_tables():
    exp, log = [0]*(2*N0), [0]*(N0 + 1)
    x = 1
    for i in range(N0):
        exp[i] = x; log[x] = i
        x <<= 1
        if x & (1 << M): x ^= PRIM
    for i in range(N0, 2*N0): exp[i] = exp[i - N0]
    return exp, log
EXP, LOG = gf_tables()

def gmul(a, b):
    return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]

def minpoly(i):
    """minimal polynomial over GF(2) of alpha^i, as an int bit mask (bit j = coeff of x^j)"""
    coset, j = [], i % N0
    while j not in coset:
        coset.append(j); j = (2*j) % N0
    poly = [1]                           # coefficients in GF(2^7), low order first
    for c in coset:
        r = EXP[c]; new = [0]*(len(poly) + 1)
        for k, p in enumerate(poly):
            new[k + 1] ^= p; new[k] ^= gmul(p, r)
        poly = new
    assert all(p in (0, 1) for p in poly)
    return sum(p << k for k, p in enumerate(poly)), frozenset(coset)

def pmul2(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r

def generator():
    g, seen = 1, set()
    for i in range(1, 2*T + 1):
        m, cs = minpoly(i)
        if cs in seen: continue
        seen.add(cs); g = pmul2(g, m)
    return g

G_POLY = generator()
DEG = G_POLY.bit_length() - 1
K0 = N0 - DEG
N, K = N0 - SHORT, K0 - SHORT
ROWS = [G_POLY << i for i in range(K)]   # shortened generator: x^i g(x), i < K, length N

def rank_masked(mask):
    rows = [r & mask for r in ROWS]
    rank = 0
    for bit in range(N):
        b = 1 << bit
        piv = next((j for j in range(rank, len(rows)) if rows[j] & b), None)
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        p = rows[rank]
        for j in range(len(rows)):
            if j != rank and rows[j] & b: rows[j] ^= p
        rank += 1
        if rank == K: break
    return rank

def p_fail(s, trials, rng):
    full = (1 << N) - 1; bad = 0
    for _ in range(trials):
        er = rng.sample(range(N), s)
        m = full
        for e in er: m &= ~(1 << e)
        if rank_masked(m) < K: bad += 1
    return bad/trials, bad

if __name__ == '__main__':
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    rng = random.Random(1)
    assert DEG == 56 and (N, K) == (120, 64), (DEG, N, K)
    out = dict(code=f'shortened BCH({N},{K}), parent ({N0},{K0}), designed d=17', trials=trials, rows=[])
    for s in list(range(12, 17)) + list(range(17, 57)):
        tr = trials if s > 16 else 2000
        pf, bad = p_fail(s, tr, rng)
        rc = min(1.0, 2.0**(-(N - K - s))) if s > 16 else 0.0
        out['rows'].append(dict(s=s, p_fail=pf, n_fail=bad, trials=tr, random_code=rc))
        print(f's={s:2d}  P_fail={pf:.3e} ({bad}/{tr})  random-code 2^-(n-k-s)={rc:.3e}', flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bch_erasure.json'), 'w'), indent=1)
