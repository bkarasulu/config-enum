# Machine-learning potentials

With `--mlip <name>`, the energies of the configurations are computed with a machine-learning interatomic potential (MLIP) through its ASE calculator. Uses include:

- ranking many configurations quickly, and passing only the best ones on to DFT;
- sequential doping to high dopant concentrations, where DFT for every configuration is too expensive;
- pre-relaxing configurations for DFT.

The model is loaded once per run. By default single-point energies are computed; `-opt` relaxes the atoms of each configuration, and `--relax_cell` the cell as well (`--fmax`, `--steps`, `--optimiser`).

```bash
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 1:3 --mlip mace --mlip_model medium-mpa-0 -opt
```

## Supported potentials

| `--mlip` | Install | `--mlip_model` (default) | Notes |
|---|---|---|---|
| `mace` (alias `mace-mp`) | `pip install mace-torch` | `small`, `medium`, `large`, `medium-mpa-0`, `medium-omat-0`, … or a path/URL to a `.model` file (`medium`) | MACE-MP foundation models or your own (fine-tuned) MACE model; `default_dtype` float64 |
| `mace-off` | `pip install mace-torch` | `small`, `medium`, `large` or a path (`medium`) | organic molecules and molecular crystals |
| `chgnet` | `pip install chgnet` | version, e.g. `0.3.0`, or a path (package default) | |
| `matgl` (alias `m3gnet`) | `pip install matgl` | a pretrained PES model name or path (`M3GNet-MP-2021.2.8-PES`) | M3GNet, TensorNet, … |
| `sevennet` (alias `7net`) | `pip install sevenn` | `7net-0`, `7net-mf-ompa`, … or a checkpoint (`7net-0`) | multi-fidelity models: `--mlip_kwargs "{'modal': 'mpa'}"` |
| `orb` | `pip install orb-models` | `orb-v3-conservative-inf-omat`, `orb-v2`, … | |
| `fairchem` (alias `uma`) | `pip install fairchem-core` | `uma-s-1`, … | needs `huggingface-cli login` and access to the UMA models; task via `--mlip_kwargs "{'task_name': 'omat'}"` (default `omat`) |
| `nequip` | `pip install nequip` | path of a compiled (nequip ≥ 0.7) or deployed (< 0.7) model (required) | |
| `allegro` | `pip install nequip-allegro` | as `nequip` | |
| `custom:pkg.module:function` | your package | passed as `model=` | `function(**mlip_kwargs)` must return an ASE calculator |

The backends are the same as in INTERFACER (the module is shared). When this version was prepared, the machine-learning path was run only through the `custom:` route with a stand-in calculator; none of the packages above could be installed in that environment. MACE is run in the `mace` job of the CI.

## Options

| Option | Meaning |
|---|---|
| `--mlip NAME` | The potential (see the table). |
| `--mlip_model MODEL` | Model name or path. |
| `--device cpu/cuda/mps/auto` | Default `auto`: the GPU if PyTorch sees one. |
| `--mlip_kwargs "{...}"` | Extra keyword arguments of the calculator, as a Python dict. |
| `-opt`, `--fmax`, `--steps`, `--optimiser` | Geometry optimisation (default: single point; fmax 0.05 eV/Å, 200 steps, BFGS). |
| `--relax_cell` | Also relax the cell (implies `-opt`). |
| `-nt N` | Run N configurations at the same time (N processes, each loading the model; mainly useful on CPUs). |

## Notes

- The structures passed on to the next step of a sequential doping are always the unrelaxed ones (dopants on the lattice sites), so that the sites and the symmetry stay well defined. The relaxed structures are written to `StepN/<formula>-config_enum-<id>.res`.
- Foundation models can be unreliable for charged or strongly off-stoichiometric defects. Check the ranking against DFT for a few configurations before relying on it, or use a fine-tuned model.
- The simple potentials (`-pot LJ/Morse/EMT`) are only for testing the set-up: LJ and Morse with ASE's default parameters do not distinguish the atom types.
