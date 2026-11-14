# Tutorial 3: the Python API

```bash
python python_api.py
```

The script shows two ways of using config_enum from Python:

1. **The building blocks**: `get_symmetry`, `get_equivalent_sites`, `select_atoms`, `enumerate_configs`, `apply_dopants`, an engine (`get_engine(...).run(atoms)`) and `config_thermo`. This gives the five symmetrically-unique Au pairs in a 32-atom Cu cell, their degeneracies (192, 48, 192, 48, 16), energies and probabilities.
2. **The workflow** (`EnumWorkflow`), which does what the `config_enum` command does and returns the results as Python objects (`results.steps[i].configs`, `results.picked`).

To use a machine-learning potential, replace the engine:

```python
engine = get_engine("mlip", CalcSettings(mlip="mace", mlip_model="small", optimise=True))
```

Any ASE calculator can be used with `AseEngine(calculator)`. See `docs/python_api.md`.
