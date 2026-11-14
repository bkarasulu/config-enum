# Tutorials

| # | Folder | System | Engine | Time | Shows |
|---|---|---|---|---|---|
| 1 | `01_Cu-Au_EMT` | Au in fcc Cu (32 atoms) | EMT (ASE) | seconds | Sequential doping, picking configurations, output files, input files |
| 2 | `02_LiCoO2_vacancies_MACE` | Li vacancies in LiCoO2 | MACE-MP-0 | minutes (CPU) | ML potentials, Wyckoff sites, several numbers of vacancies at once |
| 3 | `03_python_api` | Au in Cu | EMT (or MACE) | seconds | Python API: building blocks and the workflow |
| 4 | `04_VASP` | Li vacancies in LiCoO2 | VASP | hours (HPC) | DFT runs: INPUTS folder, Hubbard U, magnetic moments, parallel tasks, restarts |

Each folder has a README. Run a tutorial with `./run.sh` (after `pip install -e .` in the repository root) or `config_enum --input input.toml`.

Tutorial 2 needs `pip install mace-torch`.
