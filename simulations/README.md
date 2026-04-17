# FEA Simulation Suite

Self-contained C++17 verification code for the FEA architecture.
Two versions are kept here for transparency:

## `FEA_sim_v2.cpp` — current, use this

The simulation that backs the Nature Electronics manuscript. Key
properties vs v1:

- **First-principles Γ**: resonance broadening is derived from the
  lead self-energy Σ = t_c² g_L(E_F) rather than hardcoded
- **Real 2D thermal solver**: steady-state heat-diffusion SOR
  solver (not 1D slab) that includes the CMOS control-plane
  contribution
- **Contention-model crossbar**: SIM 12 uses a real queue-arbitration
  model; strided access correctly shows ~9 GOPS/zone
- **Multi-FIRE write model**: SIM 5 and SIM 10 both apply per-bit
  write fidelity from the absorption probability and automatically
  search for the minimum N_FIRE redundancy that brings per-op
  success above 99.9%. Both sims are now self-consistent
- **Langevin cross-check of Kramers**: SIM 3 now verifies that the
  first-passage-time distribution is exponential (Kramers-consistent)
- **Removed tautologies**: SIM 7 (Bernoulli self-check) and SIM 13
  (multi-hop dice-rolling) from v1 are gone; they verified only
  that the code reproduced the formula it was given

Build and run:

```bash
make v2
./FEA_sim_v2 > simulations/FEA_sim_v2_output.txt
```

or from the repo root:

```bash
make         # defaults to v2
make run     # builds and runs v2
```

## `FEA_sim_v1.cpp` — archived, do not cite

The original simulation from the initial Zenodo preprint
(DOI 10.5281/zenodo.19559255). Preserved here for transparency and
reproducibility of the preprint's reported numbers. Known
limitations compared to v2:

- Γ hardcoded at 8 meV (should be derived, ~45 meV)
- 1D thermal slab underestimates local hot spots
- Idealized crossbar throughput (no contention)
- Tautological Monte Carlo verifications in SIM 7 and SIM 13
- ALU model assumed perfect writes

Build and run:

```bash
make v1
./FEA_sim_v1 > simulations/FEA_sim_v1_output.txt
```

## Reference outputs

Both `FEA_sim_v1_output.txt` and `FEA_sim_v2_output.txt` contain
captured runs for diffing and reproducibility checks. Re-running on
a different machine should match within floating-point rounding.

## File map

```
simulations/
├── FEA_sim_v1.cpp             # archived, matches Zenodo preprint v1
├── FEA_sim_v1_output.txt      # captured v1 run
├── FEA_sim_v2.cpp             # current, backs NE manuscript
├── FEA_sim_v2_output.txt      # captured v2 run
└── README.md                  # this file
```
