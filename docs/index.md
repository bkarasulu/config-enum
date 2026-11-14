# config_enum

config_enum enumerates the symmetrically-unique configurations of point defects in a crystal (vacancies, substitutional dopants and interstitials) and ranks them by energy. It can compute the energies with:

- **VASP** or **CASTEP**, plane-wave DFT through ASE;
- **machine-learning potentials**: MACE, CHGNet, MatGL/M3GNet, SevenNet, ORB, fairchem/UMA, NequIP/Allegro, or any ASE calculator;
- **simple ASE potentials** (LJ, Morse, EMT), for quick tests.

## What it does

```text
input structure(s) ──► [interstices -ai] ──► [primitive cell -prim] ──► [supercell -r]
                                                                          │
     symmetry operations and equivalent (Wyckoff) sites (spglib) ◄────────┘      (-dry stops here)
                  │
                  ▼
     select the sites (-at, -wl, -nl, -excl) ──► all combinations of -nd dopants/vacancies (-dt)
                  │
                  ▼
     keep the symmetrically-unique configurations, with their degeneracies ──► [modifications -mod]
                  │
                  ▼
     energies (-vasp, -castep, --mlip, -pot) ──► reduced energies and probabilities at -t
                  │
                  ▼  sequential doping (-seq)
     pick the best configurations (-spick, -stype) ──► next step, starting from the picked structures
```

## Where to start

- [Installation](installation.md) and the [Quick start](quickstart.md).
- [Tutorials](tutorials.md): Au in Cu with EMT; Li vacancies in LiCoO2 with MACE and with VASP; the Python API.
- [Workflow and output files](workflow.md) and [Theory and definitions](theory.md): what is computed, and how.
- [Command-line reference](cli.md) and [Input files](input_files.md).

config_enum was originally written as a single script (`config_enum-v2.py`, kept in `legacy/`). Its options still work (except `-casm`, which was removed); see [Package layout and migration](package_layout.md) and [Fixed bugs and known issues](known_issues.md) for what changed.
