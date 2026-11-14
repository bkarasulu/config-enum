# Tutorial 6: site disorder in Li3PS4

Experimental structures often have partially occupied sites, which DFT and ML potentials cannot use directly. config_enum turns such a structure into the set of symmetrically-unique ordered models.

The Pnma structure of Li3PS4 has 16 Li sites in the cell for 12 Li atoms: `8d` (fully occupied), `4b` (occupancy 0.68) and `4c` (0.28).

```bash
pip install mace-torch
./run.sh                      # MACE-MP-0 (small), relaxing each configuration
ENGINE="-pot LJ" ./run.sh     # only to look at the configurations (seconds)
```

## What it does

1. ASE reads all sites as occupied (Li<sub>16</sub>P<sub>4</sub>S<sub>16</sub>). config_enum reports the partial occupancies:

    ```text
    Partially occupied sites in Li3PS4.cif (all sites are treated as occupied; use vacancies, -dt Va, to remove atoms):
      Li: 4 sites with occupancy 0.28 (1.1 atoms; 2.9 to remove)
      Li: 4 sites with occupancy 0.68 (2.7 atoms; 1.3 to remove)
    ```

2. Four vacancies are distributed over the eight `4b` and `4c` Li sites (`wycklist` with `aType`): 70 combinations, **16 symmetrically-unique** configurations, all Li<sub>12</sub>P<sub>4</sub>S<sub>16</sub>.

    They include the two "ordered" limits, `4b` full with `4c` empty and the reverse (degeneracy 1 each, Pnma), and the mixed occupations in between.

3. Each configuration is relaxed, and the table gives the energies, reduced energies and probabilities at 298 K. The probabilities show which occupation of the `4b`/`4c` sites is preferred, to compare with the refined occupancies.

## Things to try

- `./run.sh -t 600`: the probabilities at 600 K.
- `./run.sh -r 1 1 2 -nd 8`: the same in a doubled cell (12870 combinations).
- `./run.sh -r 1 1 2 -nd 8 -seq 2 -spick 3`: sequential doping instead of the full enumeration.
- `./run.sh --relax_cell`.

## Notes

- No reference energies are given: this tutorial could not be run with MACE when it was written. The number of configurations (16) does not depend on the potential.
- The structure is built from the data of Mercier *et al.*, Acta Cryst. B **38**, 1887 (1982).
