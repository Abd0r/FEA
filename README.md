# FEA Simulation Suite – Free Electron Absorption Architecture

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![C++17](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete simulation framework for the **Free Electron Absorption (FEA)** architecture – a radically new computing paradigm where electrons travel freely through dangling‑bond wire (DBW) pathways, and computation is performed by *Fusion Blocks* (3‑atom DB trimers) that absorb electrons on gate‑voltage command via Breit–Wigner resonance.

The simulations reproduce all quantitative claims in the accompanying paper: system clock, power density, thermal BER, Slingshot decentralised activation, 16‑atom Fusion Block ALU, and the FEA‑GP general‑purpose extension.

## 📁 Repository Structure
simulations/
├── data/ # JSON result files (generated)
│ ├── arithmetic_results.json
│ ├── ber_results.json
│ ├── btb_results.json
│ ├── halfselect_results.json
│ ├── power_breakdown.json
│ ├── sim_16atom_results.json
│ ├── slingshot_results.json
│ ├── sparse_power_results.json
│ ├── stochfp_results.json
│ ├── thermal_crosstalk_results.json
│ └── timing_results.json
│
├── fea_core.py # Fusion Core behavioral model (Python)
├── run_all.py # Master script – runs all verification modules
│
├── sim_ber.py # Kramers thermal escape + Monte Carlo BER
├── sim_timing.py # Arm/Fire/Confirm cycle timing
├── sim_power.py # Power breakdown (transit + absorption + crossbar)
├── sim_arithmetic.py # Bit‑serial FEA Word operations (AND, ADD, MUL)
├── sim_stochfp.py # Stochastic FP inference (Zone Tile Engine)
├── sim_btb.py # Branch Target Buffer (512‑entry 2‑bit saturating)
├── sim_16atom.py # 16‑atom Fusion Block ALU (AND/OR/XOR/ADD/MUL/FP16)
├── sim_slingshot.py # Slingshot hop‑chain BER and latency
├── sim_sparse_power.py # Power reduction vs. sparsity (Slingshot)
├── sim_halfselect.py # Crossbar half‑select sneak‑path Monte Carlo
├── sim_thermal_crosstalk.py # Coulomb crosstalk between adjacent DBW pathways
│
├── FEA_Chip_sim.cpp # Full chip C++ simulator (100k cores, GALS)
└── M+FEA_Chip_Sim_Output.md # Example output from the C++ simulator

text

## 🚀 Getting Started

### Prerequisites
- **Python** 3.9+ with `numpy`, `scipy` (install via `pip install numpy scipy`)
- **C++17** compiler (optional, for the full chip simulator)
- **OpenMP** (optional, for parallelisation; not required for Python scripts)

### Running the Python Simulation Suite
```bash
cd simulations
python3 run_all.py
This will execute all verification modules and produce JSON result files in simulations/data/. The console output will show PASS/FAIL for each test.

Running the C++ Full Chip Simulator
bash
g++ -std=c++17 -O2 FEA_Chip_sim.cpp -o fea_sim
./fea_sim
The C++ code simulates 100,000 Fusion Cores (6.55×10⁹ blocks) with vector ADD, matrix multiply, CAM search, sparse inference, and Slingshot messaging. The output is saved to M+FEA_Chip_Sim_Output.md.

📊 Key Verified Results
Metric	Value	Simulation
System clock	9.2 GHz	sim_timing.py
Total chip power	26.5 mW/cm²	sim_power.py
Thermal BER	8×10⁻⁷ per operation	sim_ber.py
32‑bit ADD (10⁷ parallel words)	3.49 ns, 2.87 POPS	sim_arithmetic.py
Stochastic FP (L=1024)	1.65% rel. error, 3200× speedup	sim_stochfp.py
Slingshot per‑hop latency	1.96 ns, BER < 3×10⁻¹³	sim_slingshot.py
16‑atom FP16 multiply	111 cycles, ±2 ULP	sim_16atom.py
All results are within 3% of analytical models.

🧪 Reproducibility
All random seeds are fixed (default: 42) for deterministic Monte Carlo.

JSON output files store exact parameters and results.

The paper’s analytical equations are implemented directly in the code (see fea_core.py for physical constants).

📄 Citation
If you use this code or the FEA architecture in your research, please cite:

Syed Abdur Rehman Ali. Free Electron Absorption (FEA): A New Computing Architecture of Pathways and the Fusion Block Hierarchy. Preprint, 2026.

🤝 Contributing
Issues and pull requests are welcome. Please ensure that any changes pass the full test suite (run_all.py).

📜 License
MIT License – see LICENSE file for details.

📧 Contact
Syed Abdur Rehman Ali – ra2157218@gmail.com
