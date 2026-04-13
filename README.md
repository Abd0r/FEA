# FEA — Free Electron Absorption Architecture

**A transistor-free computing architecture on hydrogen-passivated silicon.**

Electrons travel freely along dangling bond wires on H–Si(100) and are selectively absorbed by 5-atom cross-shaped dangling-bond clusters via Breit–Wigner resonance under gate-voltage control. No transistors switch in the data plane. Each Fusion Block = 1 bit. 64 Fusion Blocks = 1 64-bit Word. The CMOS control plane is an address sequencer — it never touches operand data.

> **Paper:** [`Paper/FEA-architecture.pdf`](Paper/FEA-architecture.pdf) — 17 pages, full theoretical and simulation treatment.

---

## Key Numbers

| | Value |
|---|---|
| Physical primitive | 5-atom cross DB cluster on H-Si(100) |
| 1 Fusion Block | 1 bit (absorbed or not absorbed) |
| 1 Word | 64 Fusion Blocks = 64-bit register |
| Practical block density | **3.77 × 10¹³ cm⁻²** |
| System clock (T_cycle = 108.85 ps) | **9.19 GHz** |
| Data-plane power density | 26.47 mW/cm² |
| Charging energy E_C | 0.65 eV |
| State retention τ_ret at 300 K | **52.2 ms** |
| Refresh overhead | **0.011%** (compute utilisation 99.99%) |
| ADD_64 (carry-lookahead) | **0.87 ns** (8 cycles) |
| MUL_64 (Wallace tree + CLA) | **2.18 ns** (20 cycles) |
| Slingshot hop (8-lane parallel) | **1.09 ns** |
| M4-Max-sized die (3 cm²) | **14.1 TB in-situ memory, 79.4 mW** |

**Comparison to Apple M4 Max on the same 3 cm² die area:** 110× more memory at 0.2% of the power.

---

## Simulation Suite — 13 Simulations

A self-contained C++17 simulation (`FEA_sim.cpp`, ~2,100 lines, no external dependencies) verifies the architecture end-to-end.

| # | Name | Description |
|---|------|-------------|
| 1 | 5-atom Hamiltonian + gate sweep | Jacobi diagonalisation of the 5×5 tight-binding Hamiltonian; bonding/antibonding eigenvalues vs V_NE |
| 2 | Breit–Wigner transmission | Absorption probability A(E) on/off resonance + thermal averaging |
| 3 | Kramers state retention at 300 K | τ_ret = 52.2 ms for the 5-atom cluster (E_C = 0.65 eV) |
| 4 | DBW wavepacket propagation | Crank–Nicolson on 500-site chain; v_g within 4% of analytical 2ta/ℏ |
| 5 | 64-bit ALU verification | 100,000 random trials on AND/OR/XOR/NOT/SHL/SHR/ADD/SUB/MUL — zero functional errors |
| 6 | ARM/FIRE/CONFIRM timing | 33 + 42.9 + 33 = 108.85 ps cycle |
| 7 | Thermal BER + MC verification | BER/block/cycle = 2.09 × 10⁻⁹ at E_C = 0.65 eV |
| 8 | Density, memory & power | 14.1 TB / 79.4 mW on 3 cm² die |
| 9 | Thermal map + yield MC | ΔT < 3 mK; combined yield 95% |
| 10 | General-purpose execution | Vector add, dot product, conditional branch — instruction traces |
| 11 | Room-temperature full-chip stability | M4-die epoch-based MC, 99.99% compute utilisation |
| 12 | Crossbar contention | 2,352 GOPS/zone best, 147 GOPS/zone random |
| 13 | Slingshot stress test | 64-hop chain BER 8.5 × 10⁻⁵ (below SECDED threshold) |

---

## Build & Run

Requires any C++17 compiler (clang, gcc). No external libraries.

```bash
make          # builds FEA_sim
./FEA_sim     # runs all 13 simulations
```

Or in one line:
```bash
c++ -std=c++17 -O2 -o FEA_sim FEA_sim.cpp && ./FEA_sim
```

A captured reference output from a single run is in [`FEA_sim_output.txt`](FEA_sim_output.txt) (520 lines).

---

## Architecture Overview

```
  1 Fusion Block  =  1 bit    (5-atom cross DB cluster)
  64 Fusion Blocks =  1 Word   (64-bit parallel register)
  1024 Words       =  1 Zone   (256×256 = 65,536 blocks)
  ~10⁹ Zones       =  1 Chip   (3 cm² M4-Max die)
```

**Instruction set (5 micro-ops, issued by CMOS control plane):**

- `ARM`       — address target Word within a Zone  (1 cycle)
- `FIRE`      — execute 64-bit ALU op on armed Word  (1–20 cycles)
- `CONFIRM`   — read back result via AC charge sensing  (1 cycle)
- `SLINGSHOT` — transfer 64-bit value between Words, 8-lane parallel  (10 cycles)
- `BRANCH`    — conditional jump  (1 cycle, no speculation)

---

## Key Results

**Physics.**
The 5-atom cross cluster has a smaller self-capacitance than the 16-atom 4×4 clusters used in prior work, giving E_C = 0.65 eV (vs 0.5 eV) and E_C/k_B T = 25.1 at 300 K. Kramers escape is exponentially suppressed, yielding τ_ret = 52.2 ms — approximately 330× longer retention than prior designs — which eliminates refresh overhead (<10⁻⁴).

**Arithmetic.**
Carry-lookahead adders and Wallace-tree multipliers, implemented natively on parallel Fusion Blocks, achieve 64-bit ADD in 0.87 ns and 64-bit MUL in 2.18 ns — competitive with or faster than pipelined 3 GHz CMOS.

**General-purpose execution.**
The ARM/FIRE/CONFIRM + SLINGSHOT + BRANCH instruction set is Turing-complete. Vector addition, dot product, and conditional branches are demonstrated with direct instruction traces. A compiler is future work.

**Chip scale.**
A 3 cm² M4-Max-sized die hosts 1.13 × 10¹⁴ Fusion Blocks, 14.1 TB of in-situ memory, and dissipates 79.4 mW data-plane power — roughly 504× lower power than an M4 Max GPU on the same die area. Monte Carlo simulation on the full die confirms 99.99% compute utilisation under periodic refresh at 300 K.

---

## Limitations

- **Room-temperature DB retention has not been experimentally measured.** All retention figures are Kramers-model extrapolations from cryogenic (4 K) STM measurements. This is the principal unvalidated assumption and the open experimental question.
- **No compiler.** Instruction traces in SIM 10 are hand-compiled. An LLVM backend targeting the FEA ISA is future work.
- **Fabrication.** Writing single 5-atom clusters with STM is within current capabilities; scaling to 10¹⁴ clusters on a 3 cm² die requires parallel atomic-precision patterning (directed self-assembly or template-assisted placement).

The paper describes these in detail. Experimental collaborators with atomic-precision silicon capabilities are welcome.

---

## Citation

```bibtex
@misc{ali2026fea,
  author       = {Syed Abdur Rehman Ali},
  title        = {Free Electron Absorption: A Bit-Level Transistor-Free Computing
                  Architecture on Hydrogen-Passivated Silicon},
  year         = 2026,
  howpublished = {\url{https://github.com/Abd0r/FEA}},
}
```

---

## License

Code: MIT — see [LICENSE](LICENSE).

---

**Author:** Syed Abdur Rehman Ali — Independent Researcher
[ORCID: 0009-0004-6611-2918](https://orcid.org/0009-0004-6611-2918)
