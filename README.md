# Free Electron Absorption (FEA) Architecture — Simulation Suite

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Self-contained C++17 simulation of the **Free Electron Absorption (FEA)** architecture — a transistor-free computing paradigm on hydrogen-passivated Si(100).

Electrons propagate freely through **dangling-bond wire (DBW)** pathways. Computation is performed by **16-atom (4×4) Fusion Blocks** that absorb passing electrons via Breit–Wigner resonance on gate-voltage command. No transistors exist in the data plane.

All quantitative results in the accompanying paper are reproduced by `FEA_Chip_sim.cpp` — **one file, no external dependencies**.

---

## Build & Run

```bash
git clone https://github.com/Abd0r/FEA.git
cd FEA
make
./fea_sim
```

Requires any C++17 compiler (`g++`, `clang++`). Nothing else.

---

## What the Simulation Does

### Original Analysis (`[1]`–`[10]`)

| Section | What it computes |
|---------|-----------------|
| `[1]` | Physical constants — tight-binding dispersion, Breit–Wigner, Kramers escape |
| `[2]` | 16-atom ALU: AND, OR, XOR, SHL, SHR, ADD (16-bit), MUL (16-bit), FP16 MUL |
| `[3]` | ALU statistical verification — 100,000 random trials, **zero errors** |
| `[4]` | Thermal BER — Kramers escape MC, BER vs E_C sensitivity table |
| `[5]` | Timing cycle — ARM + FIRE + CONFIRM → **9.18 GHz** system clock |
| `[6]` | Power breakdown — transit bias, absorption, gate switching → **26.47 mW cm⁻²** |
| `[7]` | Density comparison — vs 2 nm CMOS, Apple M4 GPU, NVIDIA H100 |
| `[8]` | Slingshot inter-block protocol — BER per hop, Hamming protection, 64-hop chains |
| `[9]` | Full-chip scenarios — dense, sparse (90% sparsity), neural inference, slingshot stress |
| `[10]` | Summary table of all key numbers |

### Physics Verification Suite (`SIM 1`–`SIM 6`)

| SIM | Method | Key result |
|-----|--------|-----------|
| SIM 1 | Jacobi diagonalisation of 16×16 Hamiltonian + gate sweep | 8 bonding + 8 antibonding states; B-sublattice shifts 300 meV vs 51 meV for A — independent gate control confirmed |
| SIM 2 | NEGF transmission T(E) via retarded Green's function G^r = [EI − H − Σ_L − Σ_R]⁻¹ | Multi-resonance structure confirmed; T(E_F) ≈ 0 at resonance — validates Breit–Wigner |
| SIM 3 | Crank–Nicolson wavepacket propagation on 500-site tight-binding chain | v_g = 2.24 × 10⁴ m/s — 96% of analytical 2ta/ħ |
| SIM 4 | DRAM-like refresh overhead model | 95.5% compute utilisation at τ_ret = 157.9 μs |
| SIM 5 | 32-bit multi-block carry propagation (200K trials) | Zero errors; BER consistent with single-block analysis |
| SIM 6 | Constant-interaction model, μ(N) for N = 1–16 | All Δμ/k_BT ≥ 19.3 at 300 K — 16 charge states thermally distinguishable |

### Chip-Level System Simulations (`SIM 7`–`SIM 10`)

| SIM | Method | Key result |
|-----|--------|-----------|
| SIM 7 | 2D steady-state heat diffusion (50×50 Jacobi, uniform + 10× hot-spot) | ΔT < 0.001 K; τ_ret degradation < 1%; passive cooling suffices |
| SIM 8 | Process variation Monte Carlo — 10,000 blocks, Gaussian spread in t, E_C, Γ | 92.5% parametric yield; 87.9% combined with 5% hard defects; E_C is critical parameter |
| SIM 9 | 256×256 intra-zone crossbar contention (4 access patterns) | Row-broadcast: 2,351 GOPS/zone; random: 147 GOPS/zone; 3.0 × 10⁵ TOPS chip-wide worst-case |
| SIM 10 | ARM/FIRE/CONFIRM/SLINGSHOT/BRANCH micro-instruction execution trace | Vector-add (16 elem): 1,184 cyc / 129 ns; dot-product: 5,315 cyc / 579 ns; branch: 98 cyc / 10.7 ns |

---

## Key Numbers

| Metric | Value |
|--------|-------|
| System clock | **9.18 GHz** |
| Cycle time | **108.9 ps** (ARM 33 + FIRE 42.9 + CONFIRM 33 ps) |
| Data-plane power density | **26.47 mW cm⁻²** (~3,800× lower than 2 nm CMOS) |
| Fusion Block density (practical) | **2.16 × 10¹² cm⁻²** |
| In-situ memory density | **4,325 GB cm⁻²** (volatile, 16 bits/block) |
| State retention at 300 K | **τ_ret = 157.9 μs** (E_C = 0.5 eV, Kramers model) |
| Thermal BER per operation | **6.89 × 10⁻⁷** (MC: 6.44 × 10⁻⁷, within Poisson uncertainty) |
| FP16 multiply | **111 cycles, 12.09 ns, ±2 ULP** |
| FP16 efficiency | **~10,300 TOPS W⁻¹** at 100% pathway utilisation |
| Slingshot hop latency | **1.96 ns** (18 cycles) |
| Slingshot BER per hop | **1.24 × 10⁻⁵** (Hamming-recoverable) |

---

## Architecture in Brief

```
H-Si(100) surface
│
├── Dangling Bond Wire (DBW) pathways
│     Electrons tunnel between adjacent DB atoms — no metal wires, no doping
│     Tight-binding group velocity: v_g = 2ta/ħ = 2.33 × 10⁴ m/s
│
├── 16-Atom Fusion Block  (4×4 grid, footprint 4.72 nm²)
│     A-sublattice (8 atoms, perimeter-weight) ← G_ctrl address gate
│     B-sublattice (8 atoms, interior-weight)  ← G_NE resonance gate
│     │
│     ├── State 0 (idle):    E_NE = E_F + 300 meV → T ≈ 0.9993 (transparent)
│     ├── State 1 (absorb):  E_NE = E_F           → T = 0 (quantum mirror)
│     └── Stored electron held by Coulomb blockade, E_C = 0.5 eV
│
├── ARM / FIRE / CONFIRM cycle
│     ARM    (33.0 ps)  — CMOS decoder selects zone via O(√K) crossbar
│     FIRE   (42.9 ps)  — electron injected, propagates at v_g
│     CONFIRM (33.0 ps) — AC reflectometry reads charge occupancy
│
├── O(√K) crossbar  —  256×256 intra-zone, 512 select lines per zone
│
└── Slingshot protocol  —  asynchronous DBW-native 16-bit messaging
      18 cycles/hop = 1 SOF + 16 data bits + 1 EOF
      No CMOS control-plane involvement for block-to-block data movement
```

---

## Files

```
FEA_Chip_sim.cpp    complete simulation (~2,100 lines, pure C++17)
Makefile            build: make && ./fea_sim
README.md           this file
LICENSE             MIT
```

---

## Paper

*Free Electron Absorption: A Transistor-Free Computing Architecture on Hydrogen-Passivated Silicon*  
Syed Abdur Rehman Ali — April 2026  
ORCID: [0009-0004-6611-2918](https://orcid.org/0009-0004-6611-2918)

---

## License

MIT — see [LICENSE](LICENSE).
