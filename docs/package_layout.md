# Package layout and migration guide

From version 1.0.0.dev0 the code is the installable Python package `config_enum` (in `src/config_enum`). Install it with `pip install -e .`, which also installs the `config_enum` command.

Old command lines keep working: `python config_enum.py <options>` takes the same options as the original script `config_enum-v2.py`; it calls `config_enum.cli.main`. The original script is kept unchanged in `legacy/` for reference.

## Where things moved

The functions keep their names; the main block of the script (`if __name__ == '__main__'`) became the workflow class.

| Original script | Now | Contents |
|---|---|---|
| options (`argparse`) | `config_enum/cli.py` | `build_parser()`, `main(argv)`, `parse_hubbard_u`, `parse_no_dopants`, `make_settings`, `config_from_args`, input files |
| main block: sequential loop, input-file loop, picking | `config_enum/workflow.py` | `WorkflowConfig`, `EnumWorkflow` (one method per stage), `ConfigResult`, `StepResult`, `WorkflowResults` |
| main block: atom selection, combinations, symmetry check of the configurations, doping, `-mod` | `config_enum/enumeration.py` | `select_atoms`, `parse_atom_list`, `parse_exclude_list`, `get_combinations`, `find_unique_configs`, `enumerate_configs`, `apply_dopants`, `parse_modify`, `apply_modification`, `canonical_keys` |
| main block: symmetry operations, equivalent sites; `find_interstices`, `find_equisites` | `config_enum/symmetry.py` | `get_symmetry`, `get_equivalent_sites`, `find_primitive_cell`, `site_permutations`, `invariant_symops`, `find_interstices`, `read_wyckoff_data`, `find_equisites` |
| `boltz_dist` | `config_enum/thermo.py` | `boltz_dist`, `config_thermo` |
| main block: CASM files | `config_enum/casm.py` | `prim_json`, `config_list_json`, `write_casm` |
| `handle_magmoms`, `set_hubU`, `vasp_continue`, `call_vasp_v2`, `get_potcar`, `make_potcar` | `config_enum/calculators/vasp.py` | the same, plus `prepare_vasp`, `vasp_command`, `read_pressure`, `copy_vasp_inputs` |
| `callDFT` (CASTEP part) | `config_enum/calculators/castep.py` | `call_castep`, `castep_command` |
| `callDFT` (LJ/Morse part) | `config_enum/calculators/classical.py` | `classical_calculator` |
| (new) | `config_enum/calculators/base.py`, `engines.py` | `CalcSettings`, `Engine`, `AseEngine`, `VaspEngine`, `CastepEngine`, `get_engine` |
| (new) | `config_enum/calculators/mlip.py`, `optimise.py` | `mlip_calculator` (MACE, CHGNet, ...), `relax` (shared ASE optimiser) |
| `Popen4`, `hashid`, `sorted_nicely`, `cartesian`, `combinations2`, `ncr`, `sort_dict_by_val` | `config_enum/utils.py` | the same, plus `parse_ranges` |
| `dist`, `minDist`, `volume`, `lenCell`, `lenAngleCell`, `compCells` | `config_enum/geometry.py` | the same (`lenCell_old` was removed) |
| `print` | `config_enum/log.py` | `say`, `debug` (`-v`), `warn` |
| (new) | `config_enum/io.py` | `read_struct`, `write_struct`, `write_res_result` |

Variable names of the original are kept where the code was moved: `atoms_list`, `eList`, `conf`, `degens`, `all_W`, `all_w`, `eSites`, `wyckoffs`, `nDops`, `rpt_no`, `mod_type`, `mod_value`, `mod_target`, `mod_replwith`, ...

## New command-line options

| Option | Meaning |
|---|---|
| `--mlip NAME`, `--mlip_model`, `--device`, `--mlip_kwargs` | Machine-learning potentials, see [Machine-learning potentials](mlip.md). |
| `-pot {LJ,Morse,EMT}` | The simple ASE potential used without `-vasp`, `-castep` or `--mlip` (default LJ, as before). |
| `-opt`, `--fmax`, `--steps`, `--optimiser`, `--relax_cell` | Geometry optimisation with the ASE-based engines. |
| `-ssym {current,ref}` | Symmetry used after the first step of a sequential doping, see [Theory](theory.md#sequential-doping). |
| `--wyckoff_file FILE` | Another `wyckoff.dat` for `-ai` (default: the one in the package). |
| `--input FILE`, `--save_input FILE` | [Input files](input_files.md). |
| `-q`, `--version` | Quiet output; version. |

`-od/--odir` is still accepted but has no effect (it had none in the original script either: the output goes to `StepN/input_M/`).

## Differences in behaviour

See [Fixed bugs and known issues](known_issues.md). The most important one: the original script did not recognise all symmetry-equivalent configurations, so the new version gives fewer configurations, with larger degeneracies, for the same input.
