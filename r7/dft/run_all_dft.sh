#!/usr/bin/env bash
# Reproduce every DFT result in FEA-R7 (GPAW LCAO, PBE, DZP). Long: about 3 days on an 8-core laptop.
# Results land in runs/; the paper pipeline (make -C .. all) reads them from there.
# Usage: NP=8 ./run_all_dft.sh      (needs the environment in ../environment-dft.yml)
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p runs
NP=${NP:-8}; export OMP_NUM_THREADS=1
G="mpiexec -n $NP gpaw python"
step() { echo "== $(date '+%F %T') $*"; }

step "1. slabs: 4x4 cell, 6 Si layers, one half-passivated dimer (H on left / right atom)"
python build_slab.py 4 4 6

step "2. neutral endpoint relaxation (the right endpoint is its exact mirror image)"
$G relax.py slab_4x4x6_L.extxyz runs/q0_L dzp 0 > runs/q0_L.stdout 2>&1
python mirror.py runs/q0_L.extxyz runs/q0_R.extxyz slab_4x4x6_R.extxyz

step "3. neutral CI-NEB, 3 images, IDPP start"
$G neb.py runs/q0_L.extxyz runs/q0_R.extxyz runs/neb_q0 dzp 0 3 > runs/neb_q0.stdout 2>&1

step "4. spin check at endpoint and saddle"
$G spin_fix.py > runs/spin_fix.stdout 2>&1

step "5. charged endpoints (jellium-compensated), mirrored"
for q in -1 1; do t=q$([ $q = -1 ] && echo m || echo p)
  $G relax.py runs/q0_L.extxyz runs/${t}_L dzp $q > runs/${t}_L.stdout 2>&1
  python mirror.py runs/${t}_L.extxyz runs/${t}_R.extxyz slab_4x4x6_R.extxyz
done

step "6. DB+ CI-NEB, 3 images from the neutral band"
$G neb.py runs/qp_L.extxyz runs/qp_R.extxyz runs/neb_qp dzp 1 3 runs/neb_q0_band.extxyz > runs/neb_qp.stdout 2>&1

step "7. DB- CI-NEB, 5 images, IDPP start, up to 300 steps (the 3-image band cannot follow the bridged path)"
NEB_STEPS=300 $G neb.py runs/qm_L.extxyz runs/qm_R.extxyz runs/neb_qm5 dzp -1 5 > runs/neb_qm5.stdout 2>&1

step "8. neighbour-charge tilt: extra DB on the next dimer, both bit states, q = -2 and 0"
python tilt_build.py
for q in -2 0; do for s in L R; do t=tilt_q$([ $q = -2 ] && echo m2 || echo 0)_$s
  $G relax.py tilt_$s.extxyz runs/$t dzp $q > runs/$t.stdout 2>&1
done; done

step "8b. neutral tilt in the other spin state (triplet near, broken-symmetry far); the lower state is used"
$G spin_tilt.py > runs/spin_tilt.stdout 2>&1 || echo "  note: with fixmagmom GPAW can stop after convergence while printing work functions; the triplet energy is then read from runs/spin_tilt_L.txt (see spin_tilt.py)"

step "9. configuration-coordinate single points (relaxation energy lambda)"
$G cc.py > runs/cc.stdout 2>&1
step "9b. neutral single points redone spin-polarized (1 muB start on the DB atom)"
$G cc_fix.py > runs/cc_fix.stdout 2>&1
step "done"

# ---- R7 additions (side tolerance, prefactors, bare dimers, wire) ----
step "10. harmonic TST prefactors (partial Hessian, 0.03 A steps, tight forces)"
$G htst.py q0 runs/q0_L.extxyz runs/neb_q0_band.extxyz 2 0 0.03 > runs/htst_q0_d03.stdout 2>&1
$G htst.py qm runs/qm_L.extxyz runs/neb_qm5_band.extxyz 2 -1 0.03 > runs/htst_qm_d03.stdout 2>&1

step "11. bare dimer (both H removed): symmetric and buckled starts, neutral and charged"
for q in 0 -1; do for s in sym buck; do $G bare.py $s $q > runs/bare_${s}_${q}.stdout 2>&1; done; done

step "12. bare-dimer wire: relaxed row and its pi / pi* bandwidth (hopping t = W/4)"
$G wire_band.py > runs/wire_band.stdout 2>&1

step "13. charge offset of the buckled bare dimer (Hirshfeld, site dipole)"
$G flip.py charges > runs/flip_charges.stdout 2>&1

step "14. buckling-flip barrier, charged and neutral (CI-NEB to the exact mirror image)"
$G flip.py neb -1 3 > runs/flip_neb_qm.stdout 2>&1
$G flip.py neb 0 3 > runs/flip_neb_q0.stdout 2>&1

step "15. H hop into a bare dimer from the neighbouring dimer: final-state relax, then CI-NEB"
for q in 0 -1; do for t in up down; do tag=hop_$([ $q = 0 ] && echo q0 || echo qm)_$t
  $G hop.py relax $q $t > runs/${tag}_final.stdout 2>&1
  $G hop.py neb $q $t 3 > runs/${tag}.stdout 2>&1
done; done
step "done (R7 additions)"
