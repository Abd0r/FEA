# Simulations

Three generations of the reference program, kept side by side so any claim can
be traced to the code that produced it.

```
simulations/
├── v1/                     archived reference program, v1 lineage
│   ├── FEA_sim_v1.cpp
│   └── FEA_sim_v1_output.txt
├── v2/                     archived reference program, v2 lineage
│   ├── FEA_sim_v2.cpp
│   └── FEA_sim_v2_output.txt
└── v3/                     current suite (what the paper reports)
    ├── fea_params.h        single parameter store for every v3 module
    ├── FEA_*_v3.cpp        27 modules, one per claim
    ├── outputs/            committed reference output, 28 files
    └── README.md           suite documentation and build instructions
```

## Which one is authoritative

**`v3/`.** The manuscript's numbers come from the v3 modules and are checked
against the committed files in `v3/outputs/`.

`v1/` and `v2/` are archived on purpose. They are the programs whose results
v3 reproduces, corrects or rejects, so they are kept unmodified as the
reference those comparisons are made against. Do not "fix" them: a corrected
v2 would no longer be the thing v3 is checking.

## Version history

| | What it is | Status |
|---|---|---|
| **v1** | single-file reference program | archived, read-only |
| **v2** | extended reference program | archived, read-only |
| **v3** | modular suite, one module per claim, shared parameter store | active |

Several v3 modules exist specifically to recompute a v2 figure and report
whether it reproduces. Those scenarios print both the reference value and the
recomputed one, so a divergence is visible in the output rather than hidden.

## Running it

```sh
make check          # build and run all 29 targets
make run-gamma      # one target
```

`make check` rebuilds from scratch, runs every target, and fails if any does
not pass. It writes nothing to `outputs/`, so a clean `git status` afterwards
confirms your build matches the committed reference — it does not by itself
confirm that the committed files were produced by the current sources. To
re-capture, redirect a target's own output to its file in `v3/outputs/`.

Requires a C++17 compiler. No external dependencies.
