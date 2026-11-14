# Tutorial 5: point defects in Li2S

Li2S has the antifluorite structure (Fm-3m): S on `4a`, Li on the tetrahedral `8c` sites, and empty octahedral `4b` sites. This tutorial makes three kinds of defects in the 2×2×2 supercell (96 atoms), each with one command (see `run.sh`).

```bash
pip install mace-torch
./run.sh                      # MACE-MP-0 (small), relaxing each configuration
ENGINE="-pot LJ" ./run.sh     # only to look at the configurations (seconds; the LJ energies mean nothing)
```

## (a) Two Li vacancies

```bash
config_enum -i Li2S.cif -r 2 2 2 -at Li -dt Va -nd 2 --mlip mace -opt
```

There are 2016 pairs of the 64 Li sites, and 10 symmetrically-unique ones (degeneracies 384, 192, 128, 96, 384, 128, 384, 96, 192 and 32).

## (b) Cl on a S site, with a compensating Li vacancy

```bash
config_enum -i Li2S.cif -r 2 2 2 -at S -dt Cl -nd 1 -mod "remove Li nearest 1 of Cl"
```

All 32 S sites are equivalent: one configuration (degeneracy 32). `-mod` then removes the Li atom nearest to the Cl, which gives the neutral Cl<sub>S</sub> + V<sub>Li</sub> pair (Li<sub>63</sub>S<sub>31</sub>Cl).

The eight Li neighbours of the Cl are at the same distance, and one of them is taken. They are equivalent here; when the nearest atoms are not equivalent to each other, or to study pairs at larger distances, do the two steps separately: dope first, then run config_enum on the result with `-dt Va` for the vacancy.

## (c) A Li Frenkel pair

```bash
config_enum -i Li2S.cif -ai -r 2 2 2 -at X -wl 32b -dt Li -nd 1 -mod "remove Li nearest 1 of dopant"
```

`-ai` adds the empty Wyckoff positions of the unit cell as `X` atoms: the octahedral `4b` sites and the `24d` sites between two Li atoms. With `-at X` alone, both are used:

```text
$ config_enum -i Li2S.cif -ai -at X -dt Li -nd 1 -pot LJ
 config_0     4     ...    Pm-3m (221)       (4b, octahedral)
 config_1     24    ...    Cmmm (65)         (24d)
```

In the supercell the `4b` sites are labelled `32b` (number of equivalent sites in the cell), and `-wl 32b` selects them. A Li atom is placed on one of them, and `-mod ... of dopant` removes the nearest lattice Li: a bound Frenkel pair.

## Notes

- The numbers of configurations and the degeneracies do not depend on the potential. No reference energies are given: this tutorial could not be run with MACE when it was written.
- The structure is built from the data of Buehrer *et al.*, J. Phys.: Condens. Matter **3**, 1055 (1991) (a = 5.76 Å).
