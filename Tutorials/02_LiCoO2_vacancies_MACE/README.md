# Tutorial 2: Li vacancies in LiCoO2 with MACE

This tutorial enumerates one, two and three Li vacancies in a 2×2×1 supercell of LiCoO2 (48 atoms, 12 Li) and relaxes every configuration with the MACE-MP-0 foundation model.

```bash
pip install mace-torch
./run.sh            # = config_enum --input input.toml -ow
```

The first run downloads the model. A GPU is used if PyTorch sees one (`--device`).

## What it does (`input.toml`)

1. Builds the supercell. A dry run lists the sites:

    ```text
    $ config_enum -i LiCoO2.cif -r 2 2 1 -dry -dt Va -nd 1
    Space group= R-3m (166) with tol=1.00e-03
    12a: ['Co0', ..., 'Co11']
    12b: ['Li12', ..., 'Li23']
    24c: ['O24', ..., 'O47']
    ```

2. Selects the Li atoms on the `12b` sites (`aType`, `wycklist`).
3. Enumerates the vacancies for `noDopants = "1:3"` in one step (no sequential doping):

    | Vacancies | Combinations | Symmetrically-unique | Degeneracies |
    |---|---|---|---|
    | 1 | 12 | 1 | 12 |
    | 2 | 66 | 3 | 36, 18, 12 |
    | 3 | 220 | 6 | 24, 72, 72, 36, 12, 4 |

4. Relaxes the atoms of the ten configurations with MACE (`optimise`) and prints the table with the energies, reduced energies, probabilities and space groups. The probabilities are normalised for each number of vacancies separately.

The numbers of configurations and the degeneracies above do not depend on the potential. The energies depend on the model, so no reference values are given here; this tutorial could not be run with MACE when this version was prepared (it is run in the `mace` job of the CI).

## Things to try

- `./run.sh --mlip_model medium-mpa-0`: a larger model.
- `./run.sh --relax_cell`: relax the cell as well.
- `./run.sh -nd 6 -seq 1 -spick 3`: six vacancies by sequential doping.
- `./run.sh -mod "replace Co nearest 1 of dopant with Ni"`: change the Co atom nearest to each vacancy.
- `./run.sh --mlip chgnet` (after `pip install chgnet`): another potential.
