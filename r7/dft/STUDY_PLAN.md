# Path 1: Hydrogen hop on a half-passivated dimer in H-Si(100)-2x1, by charge state and field

## Question

What is the barrier for a single hydrogen atom to hop between the two atoms of one dimer on an otherwise H-terminated Si(100)-2x1 surface? We want it for the dangling bond neutral (DB0), negative (DB-) and positive (DB+). We also want to know how much a perpendicular field, or a nearby charge, makes the two configurations unequal, since that tilt is what would let a write choose 0 or 1.

Measured anchor: H on clean Si(001), with an intra-dimer hop of 1.01 ± 0.05 eV and prefactor 10^(10.3 ± 0.5) Hz (Hill, Freelon & Ganz, PRB 60, 15896 (1999); to be confirmed from the PDF). No measurement exists for the H-terminated environment or for charged DBs.

## Model

- **Slab:** Si(100)-2x1:H with a 4x4 surface cell (8 dimers in 2 rows) and 6 Si layers. The bottom 2 layers and the bottom H atoms are fixed, with about 15 Å of vacuum.
- **Defect:** remove one top H, leaving one DB on one atom of a dimer.
- **Charge states:** 0, -1 and +1. The compensating charge is a jellium slab in the vacuum (GPAW `JelliumSlab`), with a dipole correction.
- **Engine:** GPAW in LCAO mode, PBE, Gamma point.
  - Screening uses the SZP basis on 4 layers.
  - Final numbers use the DZP basis on 6 layers, plus a 2x2 k-point check at the endpoints and saddle.
- **Path:** climbing-image NEB with 5 images, relaxed to 0.05 eV/Å.
- **Tilt:** a constant perpendicular field of ±0.1 to 0.5 V/Å (GPAW `ConstantElectricField`), and a second DB placed 0.77 to 1.5 nm away in either charge state. We report the energy difference between the two configurations and the forward and backward barriers.
- **Prefactor (optional):** harmonic transition-state theory from finite-difference vibrations of the H atom and its two Si atoms, at the minimum and at the saddle.

## Validation

1. **Benchmark against the measurement:** the H-hop barrier on the clean surface (Hill et al.), in a c(4x2) cell, compared with 1.01 eV. It is run only if time allows; otherwise we compare with published DFT values.
2. **Convergence:** layers (4 vs 6), basis (SZP vs DZP) and k-points (Gamma vs 2x2) at the endpoints and saddle. The target is barrier changes below 0.05 eV.

## Outputs

- Barrier table: E_a for DB0, DB- and DB+, plus the tilt and the asymmetric barriers against field and against a neighbour charge.
- Lifetime at 300 K for each case, using the measured prefactor range.
- Paper: "Charge-state and field control of hydrogen hopping on half-passivated Si(100) dimers", aimed at PRB or J. Phys. Chem. C.

## Compute budget

The machine is an Apple M2 with 8 cores and 16 GB RAM. An estimated 150 to 160 atoms in LCAO-DZP at the Gamma point is roughly 10 to 20 minutes per SCF; that estimate is to be measured in the validation step. Each NEB is about 30 to 60 SCF rounds of 5 images, so a DZP NEB takes about a day. Screening with SZP on 4 layers is several times faster. All runs go in the background.
