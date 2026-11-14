# Changelog

## 1.0.0.dev0 (refactor)

**Package**
- Installable `config_enum` package with a `config_enum` command. The script `config_enum.py` takes the same options as the original `config_enum-v2.py`, which is kept in `legacy/`.
- The main block of the script is split into the stages of `config_enum.workflow.EnumWorkflow`; the functions keep their names and are grouped in modules (see [Package layout](package_layout.md)).
- Results are also written to `results.json`.
- Common engine interface for VASP, CASTEP and ASE calculators.

**New features**
- Machine-learning potentials (`--mlip`): MACE, MACE-OFF, CHGNet, MatGL/M3GNet, SevenNet, ORB, fairchem/UMA, NequIP, Allegro, or a custom ASE calculator; geometry optimisation with `-opt`, `--relax_cell`, `--fmax`, `--steps`, `--optimiser`.
- Input files (`--input` TOML/YAML/JSON) and `--save_input`.
- Several numbers of dopants in one run (`-nd 0:3`), with separate statistics.
- `-mod`: several statements separated by `;`, and `dopant` as the reference also for vacancies.
- Sequential doping with `-nl`, with several input files and with several picked configurations; equivalent structures are merged within a step and not picked twice (`-ssym`).
- Interstitial sites (`-ai`): the Wyckoff positions are hard coded (`config_enum/wyckoff_positions.py`), no `wyckoff.dat` file is needed; works with `-wl` and `-seq`; only parameter-free Wyckoff positions are added.
- Partially occupied sites of the input file are reported.
- New options: `-pot`, `-ssym`, `--wyckoff_file`, `-q`, `--version`.

**Bug fixes** (see [Fixed bugs and known issues](known_issues.md))
- Symmetry-equivalent configurations were not all recognised (atom-index offset).
- Symmetry used in the later steps of a sequential doping.
- Configurational thermodynamics (free energy, entropy limits).
- VASP k-points overwritten by a Γ-point grid.
- `-excl`, `-ase`, `-stype R`, `-nd` ranges, POTCAR creation, and crashes with current ASE/spglib.

**Performance**
- The symmetry reduction works on site permutations instead of comparing coordinates of every pair of configurations: 560 triples in a 16-site cell take 0.01 s instead of 10 s.

**Removed**
- CASM support (`-casm`).

**Other**
- `p_tqdm`/`tqdm` are no longer needed.
- Added documentation, tutorials, the test suite and CI.
