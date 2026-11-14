# Installation

config_enum needs Python ≥ 3.10. Its core dependencies are NumPy, SciPy, ASE (≥ 3.23) and spglib.

```bash
git clone https://github.com/bkarasulu/config-enum.git
cd config-enum
python -m venv .venv && source .venv/bin/activate     # or a conda environment
pip install -e ".[test]"
python -m pytest                                       # ~10 s; no DFT code needed
config_enum --version
```

The `-e` (editable) install keeps the package linked to the source tree. The `config_enum` command is installed, and `python config_enum.py` works too, even without installing.

The original script needed `p_tqdm` and `tqdm`; they are no longer used.

## Optional extras

| Extra | Installs | For |
|---|---|---|
| `mace` | mace-torch | MACE-MP / MACE-OFF foundation models, your own MACE models |
| `chgnet` | chgnet | CHGNet |
| `matgl` | matgl | M3GNet, TensorNet, … |
| `sevennet` | sevenn | SevenNet |
| `orb` | orb-models | ORB |
| `fairchem` | fairchem-core | fairchem / UMA (also needs `huggingface-cli login`) |
| `nequip`, `allegro` | nequip, nequip-allegro | NequIP / Allegro models |
| `yaml` | pyyaml | YAML input files (TOML and JSON work without it) |
| `docs` | mkdocs-material, mkdocstrings | building this documentation (`mkdocs serve`) |
| `test`, `lint` | pytest, pyyaml / ruff | development |

Example: `pip install -e ".[mace,yaml]"`.

!!! note "PyTorch and GPUs"
    The ML packages pull in PyTorch. For GPU runs, install the PyTorch build that matches your CUDA version first; see [pytorch.org](https://pytorch.org/get-started/locally/). `--device auto` (the default) uses the GPU when PyTorch sees one.

## DFT codes

config_enum drives VASP and CASTEP through ASE's calculators, so they have to be set up as ASE expects:

**VASP**

- `VASP_PP_PATH` must point to the directory that contains `potpaw_PBE/`.
- `INCAR` and `KPOINTS` go in the input folder (`-id`, default `INPUTS/`).
- Give the parallel options `-np/-nt/-nn/-mpi/-vexe`, which build e.g. `srun -n 32 -N 1 vasp_std`, or set `VASP_COMMAND`.

**CASTEP**

- Set `CASTEP_COMMAND` (e.g. `mpirun -n 32 castep.mpi`), or use `-np/-mpi/-cexe`, and `CASTEP_PP_PATH`.
- Parameters can come from `<seed>.param` and `<seed>.cell` files with `-iseed <seed>`.

See [DFT codes](dft_codes.md) for details.

## Interstitial sites

The `-ai` option needs a data file with the Wyckoff positions of the space groups (`wyckoff.dat`), as the original script did. It is not distributed with the code. Put it in your home directory, or give its location with `--wyckoff_file` or the `CONFIG_ENUM_WYCKOFF` environment variable. The format is described in [Workflow](workflow.md#interstitial-sites).
