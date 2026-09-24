# FEA V3 Simulation Design

## Why this design exists

V2 asserted its headline numbers in prose. Reviewers recomputed four of them and got different
answers. V3 must therefore invert the structure: no headline number is written down, every headline
number is computed from a stated parameter set, and a unit checker fails the run if a derived value
carries the wrong dimension or the wrong order of magnitude.

This document maps each reviewer objection to one module, one equation, one output, and one
pass or fail gate. A reviewer comment with no module is an open gap, not a silent omission.

## Design rules

| Rule | Meaning |
|---|---|
| Compute, do not assert | No literal headline value in source. Derive from parameters. |
| One parameter store | One file owns every physical constant and architectural parameter. Modules import it. |
| Units on every field | Each parameter and derived quantity declares its unit. The checker converts and compares. |
| Falsifiable gate | Each module states a condition that fails the run. |
| Label the epistemic status | Every output is tagged `derived`, `estimated`, `proposed`, or `open`. |
| Uncertainty is an output | Sensitivity sweeps are first-class results, not footnotes. |

Epistemic labels, reused in the manuscript:

| Label | Meaning |
|---|---|
| derived | produced by a versioned model from stated parameters |
| estimated | calculation shown, assumptions stated, not device-validated |
| proposed | design contract exists, physical implementation unvalidated |
| open | no adequate model or measurement exists |

## Parameter store

Single source of truth: `repo/simulations/fea_params.h`.

```text
device: a, t, tc, Ec, hbar, kB, e, T, nu0, Gamma_lead, Gamma_two
geometry: die_area_cm2, zone_side_blocks, block_pitch_nm, block_bits
control: vdd, decoder_transistors, decoder_area_um2, decoder_event_J,
         pll_group_K, pll_power_W, sequencer_transistors, sense_transistors,
         sense_zone_W
timing: t_arm_ps, t_fire_ps, t_confirm_ps, vg_mps, seg_len_um, vsig_frac
io: v_bias_V, i_pathway_A, n_path
```

No module may declare its own copy of a constant. If two modules disagree, the build fails.

## Module map

Each row is one simulator. `Gate` is the condition that fails the run.

| # | Module | Reviewer objection | Computed quantity | Gate |
|---|---|---|---|---|
| M1 | `fea_budget_v3` | R5 #1, R6 #1, R7 #1: PLL 3.3 W vs 3300 W; 15 fJ gives 0.14 uW not 137.85 uW | Per-Zone and chip totals for decoder, PLL, sequencer, sensing, data plane, I-O, PDN | Order-of-magnitude unit check on each term. Any term off by 1000x fails. |
| M2 | `fea_floorplan_v3` | R6 #2, R7 #1: 12 um² decoder vs 0.174 um² per Zone; 1.7e9 x 12 um² = 204 cm² | Usable area after decoder, sensing, PLL, interconnect, PDN, FZC, edge pathways | Audits V2's own claims: `require(total > die)` must HOLD, because demonstrating that per-Zone CMOS overcommits a 3 cm^2 die IS the finding (204 cm^2 against 3.000 cm^2, and usable footprint driven negative). Whether OUR design fits its die is gated in M15, not here. |
| M3 | `fea_gamma_v3` | R5 #3, R6 #4, R7 #2: 45 meV derived vs 8 meV used | One Gamma from lead self-energy, then absorption and capture at that Gamma | Both reported Gamma values must come from one derivation. Divergence fails. |
| M4 | `fea_retention_v3` | R4 #1, R5 #4, R6 #6, R7 #4: Kramers extrapolation, 2.1 ms vs 5.31 ms at 330 K | Kramers table over Ec, nu0, T with the reviewer's own parameters reproduced | Must reproduce the reviewer's 5.31 ms at 330 K from stated inputs. Mismatch fails. |
| M5 | `fea_crosstalk_v3` | R7 #3: blocks 1.15 nm apart, neighbour Coulomb and tunnelling | Write selectivity, retention, multi-FIRE margin under occupied neighbours | Valid/invalid margin must stay positive. Overlap fails. |
| M6 | `fea_restoration_v3` | R4 #7: no gain in a transistor-free data plane, SNR collapse over 3 cm² | Link attenuation, restoration endpoint spacing, cascade output margin | Cascaded stage must clear next-stage input. Failing cascade fails. |
| M7 | `fea_multifire_v3` | R4 #10, R5 #5, R6 #7, R7 #5: P_abs 0.46 vs ideal 1.0; independent Bernoulli assumed | Single-pass capture distribution, then correlated-error sweep, then program pass rate | The correlation ceiling must bite: fidelity under common-mode blocking must fall below the independent model, and `fires_for_target(0.9999)` must become unreachable. Whether the hand-compiled program rates reconcile with this capture model is REPORTED, not gated -- `multifire` prints `program rates unreconciled` as a finding, and that divergence is carried as an open item rather than a pass condition. |
| M8 | `fea_secded_v3` | R6 #7: SECDED asserted but never simulated | Area, power, latency, refresh, correction, residual failure for the chosen code | The code parameters, the1-bit correction claim and the multi-bit residual failure rate are gated inside M8. **OPEN: the correction overhead is NOT yet folded into the M1 power or M2 area totals** -- no SECDED term appears in `budget.txt` or `floorplan.txt`. This row records that gap instead of asserting a gate that does not exist. |
| M9 | `fea_clock_v3` | R4 #10: 9.19 GHz needs justification | Latency assembled from transport, actuation, sensing, restoration, arbitration, controller | f_sys must be recomputed from summed phases. Hardcoded 9.19 fails. |
| M10 | `fea_bandwidth_v3` | R6 #5: 133 GOPS/Zone already equals chip bandwidth | Per-Zone op rate to chip aggregate through shared routing, with saturation | Dimensional check: chip bandwidth must not equal one Zone's rate by accident. |
| M11 | `fea_compare_v3` | R4 #3, #5, #6; R3 #2: M4 Max not apples-to-apples, need RC, need normalized | Energy per bit stored, area per bit, energy per confirmed operation, energy per transported bit, against DRAM, HBM, CIM, PIM, RC | Comparison must share one boundary definition. Mixed boundaries fail. |
| M12 | `FEA_fabrication_v3` (`make run-fabrication`) | R4 #8; R3 #1; CF14 #6: 1e14 clusters, STM throughput, years to pattern | Patterned atoms per second, cluster count, wall-clock time, defect and yield projection | Fabrication time must be printed as a real number, not called an engineering challenge. |
| M13 | `fea_refresh_v3` | R6 #6, R7 #4: retention sets refresh overhead | Refresh period, refresh bandwidth, refresh power, occupancy of Slingshot by refresh | Refresh IS counted inside M1 (`budget.txt` carries a refresh term) and its traffic is asserted inside the M10 ceiling. The claim that it is also counted inside M9 is **not implemented**: no clock-side refresh term exists, so that half of the original criterion is withdrawn rather than gated. |
| M14 | `fea_fzc_e2e_v3` | Control-plane replacement, already reviewed | Full transaction path with computed gates | Every gate computed. Boolean skip fails. |
| M15 | `fea_layout_v3` | No module took die area as an input; every result assumed 3 cm^2. Added to expose how layout scales | Zone grid, FZC footprint, boundary ring, clock feeds and worst-case hops as functions of die area, and capacity solved under `cells x routing + support = die` | Data-Block area must match V2's independently stated Block count within 1%. Reserving zero support must fit the die; reserving every declared structure must not. Routing slack is algebraically `sqrt(routing)`, so it is printed as an identity and never gated. |
| M16 | `fea_fzc_program_v3` | Claim owned by DECISIONS 2026-09-22, "A Zone is programmable, not fixed-function". Answers R4 #9 and R3 #2, which ask FEA to be distinguished from PIM, CIM and RC | One 0.5 cm^2 array driven by three programs: CPU-like dependency chain, GPU-like independent chains, NPU-like systolic tiling. Reports ops, cycles, active Zones, hops, ops/cycle and independent instances per die | The hardware signature must be identical across all three roles AND the three runtimes must differ. If a parameter differs between roles the result came from silicon; if the runtimes match the program was never dispatched. |
| M17 | `fea_recovery_v3` | No reviewer objection. This closes our OWN OPEN item from FZC-v0 peer recovery, where quorum latency, recovery traffic and recovery-vs-decay were all unmeasured | Seed size from the FZC ledger, the recovery latency budget against retention tau, quorum false-positive and detection probability under swept congestion, heartbeat bandwidth across every Zone, and invariant 12 epoch fencing | `T_detect + T_quorum + T_reconstruct` must be strictly below tau AND the refresh interval itself must FAIL the gate, otherwise detection is not the binding constraint. Quorum 4-of-5 must strictly reduce false recoveries against 3-of-5. The gate must move less than 5% across six decades of the unsourced latency. Heartbeat traffic must stay under 5% of refresh traffic. |
| M18 | `fea_fzc_rescue_v3` | No reviewer objection. Closes the AREA half of our own OPEN item: FZC-v0 specifies a rescue port but nothing prices the path that must reach `2.86e8` Zones, and DESIGN-V3's gate `no per-Zone CMOS scaling` demands a floorplan of where all remaining CMOS lives | Three candidate mechanisms sized against the routing budget V2's own 2x overhead leaves: a dedicated boundary tree (scales with the grid), a per-Zone CMOS decoder (scales with Zone count), and a tag detector riding the neighbour Slingshot link that already exists | The documented FZC ledger split must still sum to 535. A per-Zone decoder must EXCEED the die, or the reviewers' objection does not apply to V3 either. The boundary tree must be viable somewhere in the pitch sweep. A CMOS tag detector must NOT fit its declared width at a 5% routing budget, and the leftover must be a small positive fraction so the structural objection stands without that threshold. An in-fabric receiver must add zero area, because its Blocks are already charged. |
| M19 | `fea_thermal_v3` | Reviewer 7 #1: *revise the downstream performance AND THERMAL claims accordingly* after the control-plane budget was corrected; R4 #4 asks what the control plane costs at system level. V3 had no thermal result at all until now | 2D steady-state sheet-conduction solve over the 0.5 cm^2 die with adiabatic edges and an ideal back-face sink, heat split by LOCATION (array terms spread, boundary terms in the perimeter ring), plus d(tau)/dT derived from Kramers and three sensitivity sweeps | Solver must reproduce the textbook 1D result `t*P/(A*k)` for a uniform source to within 0.1%, or every later number is worthless. `array_power_W()` must reproduce M1's printed floor `0.023384 W` to within `1e-5` absolute (M19 currently returns `0.023390`, a difference of `6e-6`), or the two modules disagree about what the floor is. In-fabric terms alone must keep retention loss under 10%. The hot spot must land on the boundary, not the centre. The per-Zone control plane must cost over 90% of retention and over 50 K. The sweeps must contain BOTH a passing and a failing case, or the 10% criterion is not biting. |

## Required cross-module couplings

Some objections cannot be closed by one module. These couplings are mandatory:

| Coupling | Reason |
|---|---|
| M1 <- M8 | SECDED power and area must enter the chip total |
| M1 <- M13 | Refresh power must enter the chip total |
| M2 <- M8, M14 | SECDED and FZC blocks must enter the floorplan |
| M4 -> M13 | Retention sets refresh period |
| M4 -> M7 | Retention during a multi-FIRE window sets correlated-error risk |
| M6 -> M9 | Restoration endpoint spacing sets controller latency |
| M9 -> M10 | Cycle time sets aggregate bandwidth |
| M16 <- M9, M13 | M16 reports cycles only; cycle time and refresh duty stay owned upstream |
| M5 -> M7 | Neighbour occupation perturbs single-pass capture |
| M11 <- M1, M2 | Normalized comparison needs final energy and area after all overhead |
| M17 <- M4, M10 | Retention tau sets the recovery gate, and the fabric ceiling sets reconstruction time |
| M18 <- M1, M15 | The grid and die area come from M15, and the per-Zone decoder site is the same declared CMOS M1 audits |
| M19 <- M1, M4 | Heat is the power terms M1 computes, and retention against which it is judged is M4's |

## Reviewer-to-module coverage check

A reviewer comment is covered only when its module ran and its gate passed.

| Comment source | Objection | Module |
|---|---|---|
| R4 #1 | Room-temperature retention needs theory or kMC | M4, M13 |
| R4 #2 | Literature review must situate FEA vs beyond-CMOS | literature, not simulated |
| R4 #3 | M4 Max comparison conflates systems | M11 |
| R4 #4 | CMOS control-plane power overhead | M1 |
| R4 #5 | 3.8 W omits I-O energy | M1 |
| R4 #6 | Normalized core vs peripheral comparison | M11 |
| R4 #7 | Signal restoration and gain | M6 |
| R4 #8 | Manufacturing throughput | M12 |
| R4 #9 | Distinguish Fusion Memory from PIM/CIM | M11 |
| R4 #10 | Justify 9.19 GHz | M9 |
| R5 #1 | PLL arithmetic | M1 |
| R5 #3 | Gamma inconsistency | M3 |
| R5 #4 | Retention extrapolation | M4 |
| R5 #5 | Write fidelity justification | M7 |
| R5 #6 | 14.1 TB scaling justification | M2 |
| R3 #1 | Fabrication steps and cost | M12 |
| R3 #2 | RC comparison | M11 |
| R6 #1 | Power unit conversion errors | M1 |
| R6 #2 | Density omits CMOS peripherals | M2 |
| R6 #3 | Irreversible capture mechanism | open-system chain, FZC spec |
| R6 #4 | Gamma inconsistency | M3 |
| R6 #5 | Bandwidth aggregation | M10 |
| R6 #6 | Retention sensitivity and refresh | M4, M13 |
| R6 #7 | Program success vs SECDED not simulated | M7, M8 |
| R7 #1 | Control-plane power and area budget | M1, M2 |
| R7 #2 | Gamma consistency | M3 |
| R7 #3 | Neighbouring block interaction | M5 |
| R7 #4 | Room-temperature evidence and Kramers recheck | M4 |
| R7 #5 | Multi-FIRE independence assumption | M7 |
| CF14 #1, #2 | PLL power | M1 |
| CF14 #3 | Gamma | M3 |
| CF14 #4 | Retention | M4 |
| CF14 #5 | Write fidelity | M7 |
| CF14 #6 | 14.1 TB scaling | M2 |
| review.pdf #1 | Room-temperature retention extrapolated, and Haider framing | M4, M13, `research/fzc-physical-mechanism-report.md` |
| review.pdf #2 | Gap between near-unity absorption and `P_abs = 0.46` | M7 |
| review.pdf #3 | All validation is self-referential, no independent replication | open, requires a third party |
| review.pdf #4 | M4 Max comparison should be reframed or dropped | M11 |
| review.pdf #5 | `10^14` cluster manufacturing asserted not argued | M12 |
| review.pdf #6 | PLL and decoder estimates are unreferenced | M1 |

Known gaps with no simulator, carried as open rather than hidden:

| Gap | Why it is open |
|---|---|
| Ab initio or kMC stability of the 5-atom cluster at 300 K | Needs electronic-structure code and a phonon spectrum we do not have |
| Experimental room-temperature charge retention | No measurement exists |
| Independent replication of the simulation suite | Requires a third party |
| Physical Slingshot link with real attenuation | No device model |
| CMOS PDK numbers for 2 nm | No public PDK |
| Physical rescue port for peer recovery | Specified in FZC-v0 (stateless, boundary-powered, four primitives) but no device implements it. M18 resolves where it is PAID FOR: it must be in-fabric, because a CMOS detector costs 91.6% of the routing budget. The mechanism, how a stateless Block detects a rescue tag, stays open |

## Build order

| Phase | Modules | Why first |
|---|---|---|
| 1 | M1, M2 | The arithmetic errors are fatal and cheapest to fix |
| 2 | M3, M4 | Both are internal consistency defects reviewers found by hand |
| 3 | M7, M8, M13 | Reliability claims the reviewers say are unsimulated |
| 4 | M5, M6 | Physics objections needing coupled models |
| 5 | M9, M10 | Derived from phases and bandwidth once the above exist |
| 6 | M11, M12 | Comparison and fabrication depend on final energy and area |
| 7 | M14, M15, M16, M17, M18, M19 | Control-plane replacement, layout, programmability, recovery, rescue path and thermal, all added after the original 14 modules |

## Verification contract

For each module:

1. Build with `-std=c++17 -O2 -Wall`, no warnings.
2. Run, must print `PASS`.
3. Unit checker must confirm every derived quantity's dimension.
4. A reviewer-reproducible value must be reproduced. Where a reviewer recomputed a number
   by hand, the module must reproduce that hand calculation, then show the corrected value.
5. Captured output committed under `repo/simulations/outputs/`.
6. Evidence row added to `EVIDENCE.md`.

No module may report a chip-level power, area, density, bandwidth, or clock number until its
upstream couplings have run. The dependency graph above enforces this.
