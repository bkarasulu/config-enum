# config_enum

config_enum enumerates the symmetrically-unique configurations of point defects (vacancies, substitutional dopants and interstitials) in a crystal, and ranks them by energy. The energies come from VASP or CASTEP, from machine-learning interatomic potentials (MLIPs) such as MACE, or from any ASE calculator.

It runs the following workflow:

1. Read the structure; optionally add the interstitial sites (`-ai`), reduce it to the primitive cell (`-prim`) or build a supercell (`-r`).
2. Find the symmetry operations and the equivalent (Wyckoff) sites with spglib (`-dry` lists them).
3. Select the sites to dope by atom type (`-at`), Wyckoff site (`-wl`) or atom index (`-nl`).
4. Enumerate all ways of placing the dopants/vacancies (`-dt`, `-nd`) on these sites, and keep one configuration of each set of symmetry-equivalent ones, with its degeneracy.
5. Optionally modify each configuration (`-mod`), e.g. to remove the nearest cation for charge balance.
6. Compute the energies (`-vasp`, `-castep`, `--mlip`, or a simple ASE potential) and the configurational thermodynamics (reduced energies and probabilities at temperature `-t`).
7. With sequential doping (`-seq`), pick the best configurations (`-spick`, `-stype`) and repeat from step 2 with them until the target number of dopants is reached.

It can also write the input files of CASM (`-casm`).

## Installation

```bash
git clone https://github.com/bkarasulu/config-enum.git
cd config-enum
pip install -e ".[test]"      # installs the `config_enum` command and the Python package
python -m pytest              # optional: run the test suite (~10 s, no DFT code needed)
```

Machine-learning potentials are optional extras. Install the ones you need, e.g. `pip install -e ".[mace]"`. The available extras are `mace`, `chgnet`, `matgl`, `sevennet`, `orb`, `fairchem`, `nequip` and `allegro`.

To run DFT, set up the code as ASE expects:

- **VASP:** set `VASP_PP_PATH`, put `INCAR` and `KPOINTS` in the `INPUTS/` folder (`-id`), and either pass the parallel options (`-np`, `-nt`, `-nn`, `-mpi`, `-vexe`) or set `VASP_COMMAND`.
- **CASTEP:** set `CASTEP_COMMAND` (or use `-np`, `-mpi`, `-cexe`) and `CASTEP_PP_PATH`.

## Quick start

```bash
# which sites are there? (dry run: prints the space group and the equivalent sites)
config_enum -i LiCoO2.cif -r 2 2 1 -dry -dt Va -nd 1

# a quick functionality test with a simple potential (seconds): two Au atoms in a 32-atom Cu cell
config_enum -i Cu.cif -r 2 2 2 -at Cu -dt Au -nd 2 -pot EMT

# H on the interstitial sites of fcc Cu (octahedral, tetrahedral, ...)
config_enum -i Cu.cif -ai -at X -dt H -nd 1 -pot EMT

# Li vacancies with a machine-learning potential (MACE-MP-0 foundation model; GPU used if available),
# relaxing each configuration
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 1:3 --mlip mace --mlip_model medium-mpa-0 -opt

# sequential doping with VASP: four vacancies, one at a time, keeping the two lowest-energy configurations at each step
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 4 -seq 1 -spick 2 -vasp -id INPUTS -np 64 -nt 2 -mpi
```

`python config_enum.py ...` accepts the same options as `config_enum`. Run `config_enum -h` for all options. All options of the original script (`config_enum-v2.py`, kept in `legacy/`) are unchanged.

Options can also go in an input file: `config_enum --input run.toml` (TOML, YAML or JSON). Command-line options override the file, and `--save_input` writes the options of a run to a file.

## Output

```text
Step1/                              one folder per step (only Step1 without -seq)
  degeneracies                      summary: degeneracy, energy, reduced energy, probability and space group of each configuration
  input_0/config_0.res, ...         the configurations (dopants on the lattice sites, not optimised), for each input structure
  input_0/config_0/                 VASP run folders (with -vasp)
  <formula>-config_enum-<id>.res    the structures after the calculation, with the energy in the title line (AIRSS style)
  Picked/                           the configurations picked for the next step
results.json                        all results in machine-readable form
```

## Tutorials

The tutorials in [`Tutorials/`](Tutorials/README.md) are:

1. Au in Cu with EMT: sequential doping (seconds).
2. Li vacancies in LiCoO2 with MACE.
3. The Python API.
4. Li vacancies in LiCoO2 with VASP (inputs and job set-up).
5. Point defects in Li2S: vacancy pairs, Cl doping with charge compensation, Frenkel pairs.
6. Ordered models of Li3PS4 from its partially occupied Li sites.
7. Li7La3Zr2O12 (LLZO) by sequential doping over 32 steps.

## Documentation

The full documentation is built with `pip install -e ".[docs]" && mkdocs serve`. The source pages are in `docs/`:

- `docs/index.md`, `docs/quickstart.md`, `docs/workflow.md`, `docs/theory.md`, `docs/cli.md`, `docs/input_files.md`, `docs/dft_codes.md`, `docs/mlip.md`, `docs/python_api.md`.
- `docs/package_layout.md`: package layout and where the parts of the original script are now.
- `docs/known_issues.md`: bugs of the original script that were fixed, changes in behaviour, and open issues.

## Licence

config_enum is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [`LICENSE`](LICENSE).

Contributions are welcome; see [`CONTRIBUTING.md`](CONTRIBUTING.md).
