# Python API

## The workflow

```python
from config_enum.calculators import CalcSettings, get_engine
from config_enum.workflow import EnumWorkflow, WorkflowConfig

engine = get_engine("mlip", CalcSettings(mlip="mace", mlip_model="small", optimise=True))
cfg = WorkflowConfig(rep=(2, 2, 1), aType=["Li"], dType="Va", noDopants=[4], seq=1, seqpick=2)
results = EnumWorkflow(engine, cfg).run("LiCoO2.cif")

for step in results.steps:
    for c in step.configs:
        print(step.step, c.file, c.degeneracy, c.energy, c.E_reduced, c.prob, c.spacegroup, c.picked)
best = results.picked            # the configurations picked in the last step
relaxed = best[0].final_atoms    # ase.Atoms after the calculation; best[0].atoms is the unrelaxed one
```

Every `WorkflowConfig` field corresponds to a command-line option; the field comments name the option. The files are written to `cfg.base_dir` (default: the current folder), as with the command.

The steps can also be called one by one: `prepare_structure`, `get_symmetry`, `select`, `generate`, `write_configs`, `compute`, `thermodynamics`, `pick`.

## Engines

All calculations go through an engine. `engine.run(atoms, name, directory)` returns `(energy, atoms, pressure)`.

```python
from ase.calculators.emt import EMT
from config_enum.calculators import AseEngine, CalcSettings, get_engine

engine = get_engine("vasp", CalcSettings(idir="INPUTS", nprocs=64, mpirun=True, hubU={"Co": {"d": 3.32}}))
engine = get_engine("mlip", CalcSettings(mlip="chgnet", optimise=True, relax_cell=True))
engine = AseEngine(EMT(), label="EMT", settings=CalcSettings(optimise=True, fmax=0.02))   # any ASE calculator
```

## Building blocks

```python
import ase.io
from config_enum.enumeration import apply_dopants, enumerate_configs, parse_modify, select_atoms
from config_enum.symmetry import get_equivalent_sites, get_symmetry
from config_enum.thermo import config_thermo

atoms = ase.io.read("LiCoO2.cif").repeat((2, 2, 1))
data = get_symmetry(atoms, tol=1e-3)
eSites, wyckoffs = get_equivalent_sites(atoms, data)            # {'12b': ['Li12', ...]}, {'12': '12b', ...}
atoms_list = select_atoms(atoms, aType=["Li"], wycklist=["12b"], wyckoffs=wyckoffs)
confs, degens = enumerate_configs(atoms, atoms_list, 2, data["rotations"], data["translations"])
mods = parse_modify("replace Co nearest 1 of dopant with Ni")
structures = [apply_dopants(atoms, conf, "Va", atoms_list=atoms_list, mods=mods)[0] for conf in confs]
thermo = config_thermo(energies=[...], T=600., omega=degens)    # P, E_reduced, E_avg, F, S, S_max
```

See `Tutorials/03_python_api/python_api.py` for a complete example, and the API pages for all functions.
