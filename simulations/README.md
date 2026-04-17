# FEA Simulation Suite

Two self-contained C++17 simulations. Both compile with a standard
compiler and have no external dependencies.

## `FEA_sim_v2.cpp`

The reference implementation. Differences from v1:

- Γ derived from the lead self-energy Σ = t_c² g_L(E_F) instead of hardcoded
- 2D steady-state heat-diffusion SOR solver (SIM 9)
- Real queue-arbitration model for the CMOS crossbar (SIM 12)
- Multi-FIRE write-fidelity model shared by SIM 5 (block level) and
  SIM 10 (Word-level VM with 100-run Monte Carlo)
- Langevin first-passage cross-check of the Kramers retention (SIM 3)
- v1 SIM 7 (Bernoulli self-check) and SIM 13 (multi-hop dice rolling)
  removed as tautologies

Build and run:

```bash
make v2          # or  make run
./FEA_sim_v2
```

## `FEA_sim_v1.cpp`

The original implementation; it matches the numbers in the Zenodo
preprint ([DOI 10.5281/zenodo.19559255](https://doi.org/10.5281/zenodo.19559255))
and is kept here for reproducibility of the preprint figures. Known
differences from v2:

- Γ hardcoded at 8 meV
- 1D thermal slab (omits 2D hot-spot behaviour)
- Idealised crossbar throughput (no contention)
- Tautological Monte Carlo in SIM 7 and SIM 13
- ALU model assumes unit write fidelity

Build and run:

```bash
make v1
./FEA_sim_v1
```

## Reference outputs

`FEA_sim_v1_output.txt` and `FEA_sim_v2_output.txt` contain captured
runs. Rerunning on a different machine should match within
floating-point rounding.

## Files

```
simulations/
├── FEA_sim_v1.cpp
├── FEA_sim_v1_output.txt
├── FEA_sim_v2.cpp
├── FEA_sim_v2_output.txt
└── README.md
```
