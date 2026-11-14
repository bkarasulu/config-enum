# Testing and contributing

```bash
pip install -e ".[test,lint,docs]"
python -m pytest            # ~10 s; no DFT code or ML package needed
ruff check .
python docs/gen_cli.py      # regenerate docs/cli.md after changing options
mkdocs build --strict
```

## The tests

| File | Covers |
|---|---|
| `tests/test_enumeration.py` | Site selection, combinations, symmetry-unique configurations (against the algorithm of the original script in `tests/legacy_reference.py`, and against ASE's structure comparison), doping, `-mod`, structure keys |
| `tests/test_symmetry.py` | Space groups, equivalent sites and Wyckoff labels, primitive cells, site permutations, interstices and `wyckoff.dat` |
| `tests/test_thermo.py` | Configurational thermodynamics, against a direct sum over all configurations |
| `tests/test_workflow.py` | The workflow with EMT/LJ: single step, several `-nd`, sequential doping, sorting, output options, CASM, reference structures, dry run, parallel tasks, failed calculations, interstitials |
| `tests/test_calculators.py` | Engines, optimiser, ML potentials (custom calculator), VASP inputs (INCAR/KPOINTS/POTCAR with placeholder POTCARs), restart helpers |
| `tests/test_cli.py` | All original options, engine selection, input files, error messages, the `config_enum.py` script |
| `tests/test_tutorials.py` | Tutorials 1 and 3; Tutorial 2 with MACE if `CONFIG_ENUM_SLOW_TESTS=1` and `mace-torch` is installed |

VASP and CASTEP are not run by the tests: only the preparation of the VASP inputs is tested. Please check the inputs of one configuration before production runs, and report problems.

## Guidelines

- Keep the existing command-line options and function names working.
- Add a test for every bug fix and new feature.
- New options go in `config_enum/cli.py` (`build_parser`, `make_settings` or `config_from_args`) and in `WorkflowConfig` or `CalcSettings`; then run `python docs/gen_cli.py`.
- A new engine is a subclass of `config_enum.calculators.base.Engine` with a `run(atoms, name, directory)` method returning `(energy, atoms, pressure)`, registered in `get_engine`.
