# FEA — Free Electron Absorption Architecture

A transistor-free computing architecture on hydrogen-passivated Si(100).
Electrons travel along dangling bond wires (DBWs); 5-atom cross-shaped
dangling-bond clusters selectively capture passing electrons via
Breit–Wigner resonance under gate-voltage control. One cluster (Fusion
Block) stores one bit; 64 Fusion Blocks form a 64-bit Word.

<p align="center">
  <img src="docs/img/fusion_block_hierarchy.png" width="90%" alt="Architectural hierarchy: 5-atom Fusion Block, 64-bit Word, and 256×256-block Zone on H-Si(100)">
</p>

<p align="center"><em>Fusion Block (a), 64-bit Word (b), 256×256-block Zone (c). A 3 cm² die contains ~1.7 × 10⁹ Zones.</em></p>

**Preprint:** [10.5281/zenodo.19559255](https://doi.org/10.5281/zenodo.19559255)
· **Paper PDF:** [`Paper/FEA-architecture.pdf`](Paper/FEA-architecture.pdf)
(22 pages, 8 figures, 14 simulations, 13 references).

The Zenodo preprint reports numbers from `FEA_sim_v1.cpp`.
The current reference implementation is
[`simulations/FEA_sim_v2.cpp`](simulations/FEA_sim_v2.cpp), which
replaces the hardcoded Γ with a first-principles derivation, adds a 2D
thermal solver, adds a contention-model crossbar, and makes the
ALU/VM simulations consistent under a multi-FIRE write-fidelity model.
Numbers below are v2.

---

## Key Numbers

| | Value |
|---|---|
| Physical primitive | 5-atom cross DB cluster on H-Si(100) |
| 1 Fusion Block | 1 bit |
| 1 Word | 64 Fusion Blocks |
| Practical block density | 3.77 × 10¹³ cm⁻² |
| System clock (T_cycle = 108.85 ps) | 9.19 GHz |
| Resonance broadening Γ (derived) | 45 meV |
| Data-plane power density | 26.47 mW/cm² |
| Charging energy E_C | 0.65 eV |
| Retention τ_ret at 300 K | 52.2 ms |
| Refresh overhead | 0.011% |
| ADD_64 single-FIRE ideal / multi-FIRE (N=18) | 0.87 ns / 2.72 ns |
| MUL_64 single-FIRE ideal / multi-FIRE | 2.18 ns / 4.03 ns |
| Slingshot hop (8-lane parallel) | 1.09 ns |
| 3 cm² die: in-situ memory | 14.1 TB |
| 3 cm² die: data-plane power | 79.4 mW |
| 3 cm² die: total chip power (incl. CMOS control plane) | ~3.8 W |

---

## Simulation Suite (v2, 14 simulations)

Source: [`simulations/FEA_sim_v2.cpp`](simulations/FEA_sim_v2.cpp)
(~1,370 lines, C++17, no external dependencies).

| # | Name | Result |
|---|------|--------|
| 1 | Hamiltonian + lead self-energy | Γ = 45 meV (derived from Σ = t_c² g_L) |
| 2 | Transmission: Green's function vs Breit–Wigner | T(E) from full GF; injection-averaged capture |
| 3 | Kramers retention + Langevin MC | τ_ret = 52.2 ms; FPT distribution exponential |
| 4 | Wavepacket with embedded cluster | On/off absorption contrast ~1,066× |
| 5 | CLA ADD with multi-FIRE redundancy | N_FIRE = 18 → 0/1,000 failures → 2.72 ns |
| 6 | ARM/FIRE/CONFIRM timing | 33 + 42.9 + 33 = 108.85 ps |
| 8 | Density, memory, power breakdown | 14.1 TB; P_transit + P_absorb + P_gate |
| 9 | 2D steady-state thermal (SOR) | ΔT_max ≈ 1.07 K with CMOS, <3 mK data-plane only |
| 10 | FEA VM — three programs, 100-run MC | vector add 96%, dot product 98%, branch 98% |
| 11 | Full-chip stability (binomial MC) | 99.99% compute utilisation under τ/2 refresh |
| 12 | Crossbar contention (real arbitration) | 147 / 133 / 9 GOPS/zone (row / random / strided) |
| 14 | Cross-die access (fat-tree routing) | 95th %ile 32 hops (34.8 ns), max 34 hops (37.0 ns) |

v1's SIM 7 and SIM 13 were tautological Monte Carlo self-checks and have been removed; the underlying physics is now covered by SIM 3 (Langevin) and SIM 5 (block-level CLA).

---

## Build & Run

C++17 compiler (clang, gcc). No external libraries.

```bash
make          # builds FEA_sim_v2
make run      # runs v2
make v1       # builds the archived v1 for reproducibility
```

One-liner:

```bash
c++ -std=c++17 -O2 -o FEA_sim_v2 simulations/FEA_sim_v2.cpp && ./FEA_sim_v2
```

Captured reference outputs:
[`simulations/FEA_sim_v2_output.txt`](simulations/FEA_sim_v2_output.txt),
[`simulations/FEA_sim_v1_output.txt`](simulations/FEA_sim_v1_output.txt).

---

## Architecture

```
  1 Fusion Block  =  1 bit    (5-atom cross DB cluster)
  64 Fusion Blocks =  1 Word   (64-bit parallel register)
  1024 Words       =  1 Zone   (256×256 = 65,536 blocks)
  ~10⁹ Zones       =  1 Chip   (3 cm² M4-Max die)
```

<p align="center">
  <img src="docs/img/cim_pim_fea.png" width="95%" alt="CIM, PIM, and FEA compared">
</p>

<p align="center"><em>Compute–memory integration. (a) CIM: array + peripheral ADCs + accumulators + separate decoder. (b) PIM: logic near DRAM banks, fetch–execute boundary preserved per bank. (c) FEA: each Fusion Block is simultaneously the memory cell and the compute unit.</em></p>

Instruction set (5 micro-ops, CMOS control plane):

- `ARM` Zone, Word — address target Word (1 cycle)
- `FIRE` Op — execute ALU op on armed Word (1–20 cycles)
- `CONFIRM` — read back result via AC charge sensing (1 cycle)
- `SLINGSHOT` src, dst — 8-lane parallel 64-bit transfer (10 cycles)
- `BRANCH` cond, offset — conditional jump (1 cycle, no speculation)

---

## Sparse-Workload Power

The data plane has no dark-silicon floor. The CMOS control plane
(~3.8 W, dominated by PLL distribution) does, unless actively
duty-cycled. The three relevant numbers at 5% activation (95%-sparse
transformer):

| Layer | Power at 5% activation |
|---|---|
| Data plane only | ~4 mW |
| Full chip, static CMOS | ~3.81 W |
| Full chip, duty-cycled CMOS (100 mW floor) | ~0.29 W |

<p align="center">
  <img src="docs/img/power_vs_sparsity.png" width="85%" alt="Full-chip power vs activation fraction">
</p>

<p align="center"><em>Full-chip power vs activation fraction (3 cm² die, log–log). Green: data plane only. Blue: total chip, static CMOS. Purple: total chip, duty-cycled CMOS. Dashed red: Apple M4 Max (40 W). Dotted orange: human brain (20 W).</em></p>

---

## Limitations

- Room-temperature DB retention has not been experimentally measured. All retention figures are Kramers-model extrapolations from cryogenic (~4 K) STM measurements.
- No compiler. Instruction traces are hand-compiled.
- Fabrication of 10¹⁴ clusters on a 3 cm² die requires massively parallel atomic-precision patterning; current STM demonstrations operate at the 10²–10⁴-atom scale.
- CMOS control-plane power is an analytical estimate; physical design would be required to confirm it and to bound clock-distribution, I/O, and PDN contributions.

---

## Citation

```bibtex
@misc{ali2026fea,
  author    = {Ali, Syed Abdur Rehman},
  title     = {Free Electron Absorption: A Bit-Level Transistor-Free Computing
               Architecture on Hydrogen-Passivated Silicon},
  year      = 2026,
  month     = apr,
  publisher = {Zenodo},
  version   = {v1.0},
  doi       = {10.5281/zenodo.19559255},
  url       = {https://doi.org/10.5281/zenodo.19559255},
}
```

Plain text:

> Ali, S. A. R. (2026). *Free Electron Absorption: A Bit-Level Transistor-Free Computing Architecture on Hydrogen-Passivated Silicon* (v1.0). Zenodo. https://doi.org/10.5281/zenodo.19559255

---

## License

Code: MIT — see [LICENSE](LICENSE).

---

Syed Abdur Rehman Ali · Independent Researcher ·
[ORCID 0009-0004-6611-2918](https://orcid.org/0009-0004-6611-2918)
