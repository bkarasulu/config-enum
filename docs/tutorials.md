# Tutorials

The tutorials are in the `Tutorials/` folder of the repository. Each has a README with the expected results and a `run.sh` script.

| # | Folder | System | Engine | Time | Shows |
|---|---|---|---|---|---|
| 1 | `01_Cu-Au_EMT` | Au in fcc Cu (32 atoms) | EMT (ASE) | seconds | Sequential doping, picking configurations, output files, input files |
| 2 | `02_LiCoO2_vacancies_MACE` | Li vacancies in LiCoO2 | MACE-MP-0 | minutes (CPU) | ML potentials, Wyckoff sites, several numbers of vacancies at once |
| 3 | `03_python_api` | Au in Cu | EMT (or MACE) | seconds | Python API: building blocks and the workflow |
| 4 | `04_VASP` | Li vacancies in LiCoO2 | VASP | hours (HPC) | DFT runs: INPUTS folder, Hubbard U, magnetic moments, parallel tasks, restarts |

Tutorials 1 and 3 are run by the test suite (`tests/test_tutorials.py`). Tutorial 2 needs `pip install mace-torch`.
