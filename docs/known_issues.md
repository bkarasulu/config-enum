# Fixed bugs and known issues

This page lists the bugs of the original script (`legacy/config_enum-v2.py`) that were fixed in the refactor, the deliberate changes in behaviour, and the issues that are still open. Regression tests are in `tests/`.

## Fixed

### Results

| Area | Problem | Fix |
|---|---|---|
| Symmetry-equivalent configurations | The check used the coordinates of the atom *before* each selected atom (`scpos[at-1]`; a comment in the script notes that it was originally `scpos[at]`). The symmetry images therefore belonged to other atoms, and many equivalent configurations were not recognised: e.g. two Au atoms on the 16 sites of a 2×2×1 Cu cell gave 18 "unique" configurations instead of 5, two Na sites in the conventional NaCl cell 2 instead of 1. Degeneracies were wrong accordingly, and with them the reduced energies and probabilities. | The positions of the selected atoms are used. The result is identical to the original algorithm with `scpos[at]` (tested), and agrees with ASE's structure comparison. |
| Sequential doping: symmetry after Step 1 | The symmetry operations of the original (undoped) input were used in all steps. Sites that are equivalent in the undoped structure but at different distances from the dopants already present were merged, so inequivalent configurations were lost. With `-prim` or `-r`, the operations were those of the input cell, not of the cell being doped. | The symmetry of the current (doped) structure is used, and structures that are equivalent through the symmetry of the Step 1 structure are merged afterwards (see [Theory](theory.md#sequential-doping)). `-ssym ref` gives the old behaviour (with the operations in the right cell). |
| Sequential doping: remaining sites | The list of sites left for the next step was updated by comparing atom positions with the last structure of the previous step; with more than one picked configuration (`-spick`) or input file, the list of one configuration was used for the others. | Each configuration carries its own list of remaining sites; `-nl`, `-mod` and several input files work with `-seq`. |
| `boltz_dist` | The upper limit of the entropy used a global variable (`degens`) instead of the argument; the free energy was computed as −k<sub>B</sub>T·ΣP instead of −k<sub>B</sub>T ln Z; the "complete space" values ignored the degeneracies. (These were only printed with `-v`; the probabilities and reduced energies were right.) | `config_thermo`, tested against a direct sum over all configurations. |
| Several numbers of dopants (`-nd 0:2`) | Crashed without `-seq` (`int('0:2')`); otherwise only the configurations of the last number were computed, with degeneracies stored at the wrong index. | All numbers are enumerated and computed in one step; the statistics is done for each number separately. `-nd 0` gives the undoped structure. |
| Excluded atoms (`-excl`) | Without `-dexcl`, the excluded atoms were still modified when `-dt Va` was used (their symbol was set to the empty `-dexcl_type`); with `-dexcl` and no `-dexcl_type` the run failed. | Excluded atoms are untouched unless `-dexcl` is given; `-dexcl_type` defaults to `-dt`. |
| `-mod ... of dopant` | The reference indices were those before the deletion of the vacancies, so the wrong atoms were used with `-dt Va`. | The positions of the doped/vacant sites are the reference. |
| `-ase` | After merging, the file list and the numbering no longer matched the configurations. | The check is done before the files are written. |
| `-stype R` | Shuffled the four columns of the table, not the configurations. | Random order of the configurations. |
| `-stype S` | (Fixed in the last original version; kept:) highest symmetry first, lowest energy first within the same space group. | |
| VASP k-points | `KPOINTS` was copied from the input folder, but ASE then wrote its own Γ-point-only file over it. | The k-points of the copied `KPOINTS` are read into the calculator. **VASP energies can differ from runs with the original script for this reason.** |
| VASP run folders | `call_vasp_v2` changed into a folder that did not exist and continued in the current one; `RUN10` was sorted before `RUN2` in `vasp_continue`. | Explicit run folders (`StepN/input_M/config_K/`); numeric order. |
| Wyckoff labels | Sets of equivalent atoms with the same label were merged and the multiplicity recomputed as (number of atoms)/(number of atom types), which gave wrong labels for two sets of the same atom type (e.g. `16c` for two `8c` sets) and with a reference structure (`-ref`) of different composition. With `-prim`, the multiplicity factor could become 0. | Each set is labelled by its own multiplicity; the factor for `-prim` is the volume ratio. The labels of the common cases (e.g. LLZO `96h`, `24d`) are unchanged. |

### Crashes

- `get_potcar` used a variable (`vasppp`) that only existed if `VASP_PP_PATH` was set when the script started.
- `find_interstices` (`-ai`): the space-group argument was overwritten while reading `wyckoff.dat`, so the option stopped at the first empty line of the file; the last block of the file was skipped; Wyckoff positions whose site symmetry is written with a space (e.g. `48 i m.m 2`) were not recognised, and their coordinates were added to the previous position.
- `-ai` with `-seq`: all unused interstices were removed after the first step, so there were no sites left for the next one.
- `-mgm`: the magnetic moments were passed to ASE as strings.
- `-hubU` with current ASE: logical INCAR tags were given as `".True."` and `1`; recent ASE needs Python booleans.
- Current ASE and spglib: `atoms.set_calculator`, `ase.spacegroup.get_spacegroup`, `ase.io.castep.read_cell` and dictionary access to the spglib dataset are deprecated or removed.
- `-nd N` with fewer than N selected sites and unreadable input files give a message instead of a traceback.

## Deliberate changes

- **CASM support was removed** (`-casm`, `prim.json`, `config_list.json`, `POS` files). In the original script the option stopped with an error after the enumeration (`outff.writelines(sbtr11)`).
- **VASP structure files** (`-ot vasp`) are written with the atoms sorted by type, so that the POSCAR has one block per element (and `-mp` one POTCAR entry per element).
- **Dependencies:** `p_tqdm` and `tqdm` are no longer needed; `-nt` uses Python's `multiprocessing`.
- **`-ref POSCAR`:** the original script set the symmetry tolerance to 1.0 when the reference file was called `POSCAR`, which makes all configurations equivalent. The tolerance is now always `-tol`.
- **`-castep` no longer forces `-stype E`.**
- **Failed calculations** are skipped (energy 0, P1) for all engines, not only VASP, and are not picked while other configurations are available.
- **Step folders:** as before, the output is always in `Step1/` (and further steps), also without `-seq`.
- **`.res` result files** now contain the energy in the title line (it was `None`).
- **Interstices (`-ai`):** only Wyckoff positions without free parameters are added. The original set the free parameters to 0, which gives arbitrary points (e.g. 96 extra sites in diamond Si from (1/8,y,-y+1/4) with y = 0). The Wyckoff positions are now part of the code (`config_enum/wyckoff_positions.py`); no data file is needed.
- **Too many combinations:** the message for more than 2×10⁴ combinations is kept as information; nothing is skipped (as before). The symmetry reduction is about three orders of magnitude faster than the original, so larger problems are feasible; all combinations are still generated in memory.

## Open issues

- **VASP and CASTEP were not run** in the environment where the refactor was done: the VASP input preparation is tested (INCAR, KPOINTS, POTCAR, magnetic moments, Hubbard U), the run and restart logic is ported from the original and untested; the CASTEP engine is ported and untested. Please test with a small system first.
- **ML potentials:** only the `custom:` route was run; see [Machine-learning potentials](mlip.md).
- **`-ai`** needs a conventional cell in the first setting of the space group in the table (origin choice 1 where there are two); structures given in another setting get interstices at the wrong places. Check the result with `-dry`-like runs on a small cell (`-pot LJ`) before production runs.
- **Wyckoff table:** the data in `config_enum/wyckoff_positions.py` are those of the `wyckoff.dat` of the original script, which was taken from VESTA; check that this is acceptable before making the repository or a release public.
- **`-ref`** must have the same cell and atom order as the structure to be doped; a warning is printed if the cells differ. With a reference of higher symmetry than the structure, configurations are merged according to the reference (by design).
- **Sequential doping after relaxation:** the next step starts from the unrelaxed picked structures. Starting from relaxed structures would need a mapping of the relaxed atoms back to the lattice sites.
- **Multiple dopant types** in one run are not supported (to-do item of the original). They can be done one after the other with `-dexcl`/`-dexcl_type`, or by running config_enum on the picked structures.
- **Charged defects / Ewald energies** (to-do item of the original) are not implemented.
- **Pre-optimisation of the input structure** (to-do item of the original) is not implemented.
- The sorted list of a step with several numbers of dopants (`-nd 0:3` without `-seq`) mixes compositions; use the per-composition tables.
