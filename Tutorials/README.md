# Tutorials

| # | Folder | System | Engine | Time | Shows |
|---|---|---|---|---|---|
| 1 | `01_Cu-Au_EMT` | Au in fcc Cu (32 atoms) | EMT (ASE) | seconds | Sequential doping, picking configurations, output files, input files |
| 2 | `02_LiCoO2_vacancies_MACE` | Li vacancies in LiCoO2 | MACE-MP-0 | minutes (CPU) | ML potentials, Wyckoff sites, several numbers of vacancies at once |
| 3 | `03_python_api` | Au in Cu | EMT (or MACE) | seconds | Python API: building blocks and the workflow |
| 4 | `04_VASP` | Li vacancies in LiCoO2 | VASP | hours (HPC) | DFT runs: INPUTS folder, Hubbard U, magnetic moments, parallel tasks, restarts |
| 5 | `05_Li2S_defects` | Li2S (antifluorite) | MACE-MP-0 | minutes | Vacancy pairs, aliovalent doping with charge compensation (`-mod`), interstitials and Frenkel pairs (`-ai`) |
| 6 | `06_Li3PS4_site_disorder` | Li3PS4 (Pnma) | MACE-MP-0 | minutes | Ordered models from partially occupied sites; `-wl` with `-at` |
| 7 | `07_LLZO_sequential` | Li7La3Zr2O12 garnet | MACE-MP-0 | ~1 h (CPU) | Sequential doping over 32 steps in the primitive cell |

Each folder has a README. Run a tutorial with `./run.sh` (after `pip install -e .` in the repository root) or `config_enum --input input.toml`.

Tutorials 2 and 5–7 need `pip install mace-torch`; with `ENGINE="-pot LJ" ./run.sh`, tutorials 5–7 run without it, to look at the configurations only (the LJ energies mean nothing).
