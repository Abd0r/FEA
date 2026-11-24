# FEA-R7: room-temperature single-electron compute-in-memory on H-Si(100)

Code and data for the manuscript **"Free Electron Absorption: room-temperature single-electron
compute-in-memory on hydrogen-terminated silicon"** (R7). Every number, figure and table in the paper
is produced by the scripts here, starting from the committed DFT results.

## Quick start (a few minutes, any laptop)

```bash
cd r7
python3 -m pip install -r requirements.txt   # numpy, matplotlib, ase
make check                                   # rebuild everything, then verify the key numbers
```

`make check` prints one PASS/FAIL line per quantity: the DFT barriers, the tilt, the capacity, the
energies, the relaxation energy λ and every thermal-margin result of the chip simulation. It exits non-zero if anything
differs. Figures go to `paper/figures/`. Every number quoted in the manuscript is written to
`paper/numbers_r7.tex` as a LaTeX macro. The manuscript source will be released here on publication.

## What runs

| Step | Script | Produces |
|---|---|---|
| Capture model (R5) | `paper/fea_r5.py` | `numbers_r5.json`, Supplementary Figs. S1–S2 |
| Line/control model (R6) | `paper/fea_r6.py` | `numbers_r6.json` |
| Chip, DFT extraction, stability, tilt | `paper/fea_r7.py` | `numbers_r7.json` |
| Retention model (5 optimisations) | `opt/retention_model.py` | Fig. 3a |
| Whole-chip simulation | `opt/fea_chip_sim.py` | `chip_sim_output.txt`, `.json` |
| Adaptive vs fixed refresh, real-world day | `opt/fea_chip_sim_adaptive.py` | `chip_sim_adaptive_output.txt`, `.json` |
| Figures | `paper/fea_r7_figs.py` | `paper/figures/r7_fig*.pdf` |
| Numbers in the text | `paper/make_macros_r7.py` | `paper/numbers_r7.tex` (LaTeX macros) |

The manuscript contains no hand-typed numbers. Each one is a macro written by `make_macros_r7.py`.

## DFT (optional, days)

All DFT inputs, relaxed structures, band images, optimiser logs and JSON summaries are in
`dft/runs/`. The full GPAW text outputs are in `dft/runs/gpaw_text_outputs.tar.xz`. To recompute
them from scratch (GPAW 26.7, PBE, DZP, about 3 days on 8 cores):

```bash
micromamba create -f environment-dft.yml && micromamba activate fea-dft
NP=8 make dft
```

`dft/run_all_dft.sh` lists every step in order: slab construction, endpoint relaxations, mirror
endpoints, CI-NEB for DB0, DB+ and DB- (5 images), the spin check, the neighbour-charge tilt in both spin states, and
the configuration-coordinate (λ) single points.

## Layout

```
r7/
  paper/   analysis scripts; figures and number macros are generated here
  opt/     retention model and chip simulations
  dft/     DFT scripts, structures and results (runs/)
  check_r7.py, reference_r7.json   reproduction check
```

Earlier revisions (v1–v3, C++) remain under `simulations/` for provenance. R7 supersedes their
numbers.
