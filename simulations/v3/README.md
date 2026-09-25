# FEA Simulation Suite

27 self-checking V3 targets plus the two legacy V1/V2 sources they audit, 29
C++17 files in total. Every target compiles with a standard compiler, has no
external dependencies, and prints `PASS` only when all of its gates hold.

Build and run one target with `make run-<target>`, for example `make run-program`.

Build and run every target with `make check` (from the repository root). It
rebuilds the whole suite, prints PASS or FAIL per target, shows the failing
output for any failure, and exits non-zero if anything did not pass — so it is
a single reproducible command.

## `FEA_fzc_e2e_v3.cpp`

One transaction computes a dual-rail score, screened bias, clocked storage,
a sensor reading, a three-replica vote, a link walk, and a commit-log check.
It commits once only if every computed gate passes. A bad codeword, weak
screening, missing clock, sensor offset, split vote, down link, or repeated
transaction id stops the path with no commit. Parameters are declared. This is
not a silicon proof.

Build and run:

```bash
make run-fzc-e2e
```

## `FEA_fzc_floorplan_v3.cpp`

A Block-count ledger for one Zone, which is 66,071 Blocks: the 65,536-Block
data array (1024 Words) plus the 535-Block FZC. It rejects a budget that omits
pathways or puts CMOS back in every Zone. Under the declared default widths,
the FZC is 535 Blocks, so the provisional 512 target fails. Cutting state to 4
bits per group fits in 512, which shows the target is conditional. Shared edge
pathways can exceed the FZC. These are allocation counts, not a drawn layout,
nanometre area, or terabyte density.

Build and run:

```bash
make run-fzc-floorplan
```

## `FEA_fzc_screened_v3.cpp`

The same four-DB geometry with a declared Yukawa screening length and a
temperature-dependent escape hold. Long screening keeps the near-rail bias.
Short screening drops it below the store margin and blocks a 3-unit cascade hop.
Low theta retains the stored occupation computed in this file. High theta empties it. `lambda`
and `E_bind` are not silicon measurements.

Build and run:

```bash
make run-fzc-screened
```

## `FEA_fzc_selector_v3.cpp`

A stated four-DB selector geometry plus a clocked rate equation. Command-rail
bias is computed from 1/r sums in declared lattice units. A sine clock is the
only energy input. Storage requires a favorable bias above `0.5`, a falling
clock, and reservoir exchange. The ideal sensor reads stored occupation; a
declared flip can change CONFIRM. Far separation fails storage. Lattice units
are not nanometres, and the rates are not a calibrated Si Hamiltonian.

Build and run:

```bash
make run-fzc-selector
```

## `FEA_fzc_opensystem_v3.cpp`

An open-system accounting gate for one FZC selector stage. It requires five
named channels (reflection, transmission, temporary actuator occupation,
stored occupation, and reservoir exchange) to conserve probability. Transmission
or a stored label without a reservoir is rejected as capture. Restoration and a
two-stage cascade are allowed only when clock/bias energy and output margin are
both declared. The numbers are accounting parameters, not device rates or energies.

Build and run:

```bash
make run-fzc-opensystem
```

## `FEA_fzc_recognition_v3.cpp`

A parameterized, dimensionless dual-rail FZC pattern-recognition model. It
measures target/invalid code separation, rejection of a corrupted rail pair,
declared neighbour-offset sensitivity, and seeded disorder-induced misses and
false fires. Its values are abstract selectivity parameters, not calibrated
DBW/Si dangling-bond device quantities. It does not model transport, capture,
relaxation, sensing, gain/restoration, timing, energy, or retention.

Build and run:

```bash
make run-fzc-recognition
```

## `FEA_fzc_reliability_v3.cpp`

The V3 protected-controller-state simulator. It uses triple replicated logical
controller state and tests single-replica correction, uncorrectable disagreement
safe halt, bounded retry, exactly-once commit, and a seeded mix of one-fault and two-fault cases.
The injected faults are deliberately not physical error-rate claims.

Build and run:

```bash
make run-fzc-reliability
```

## `FEA_slingshot_v3.cpp`

The V3 unified Slingshot network simulator. It models typed packet priority,
fixed-size ingress buffers, coordinate routing, per-link arbitration, explicit
ingress rejection, link-fault rejection, and packet-conservation invariants.
Its unit is an arbitration round, not physical time. It does not claim link
energy, bandwidth, DBW attenuation, capture, or restoration.

Build and run:

```bash
make run-slingshot
```

## `FEA_fzc_v3.cpp`

The V3 Fusion Zone Controller and unified Slingshot architectural simulator.
It models protocol semantics, not unresolved device physics. It verifies
boot, ARM, FIRE, CONFIRM, commit, retry, duplicate suppression, refresh-source
rules, recovery, and safe state transitions on a small Zone mesh. The actuator
is intentionally an abstract interface until a physical FZC-pattern-to-actuator
model exists.

Build and run:

```bash
make run-fzc
```

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
├── FEA_fzc_v3.cpp
├── FEA_slingshot_v3.cpp
├── FEA_fzc_reliability_v3.cpp
├── FEA_fzc_recognition_v3.cpp
├── FEA_fzc_opensystem_v3.cpp
├── FEA_fzc_selector_v3.cpp
├── FEA_fzc_screened_v3.cpp
├── FEA_fzc_floorplan_v3.cpp
├── FEA_fzc_e2e_v3.cpp
└── README.md
```
