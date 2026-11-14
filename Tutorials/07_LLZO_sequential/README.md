# Tutorial 7: LLZO by sequential doping

Cubic Li7La3Zr2O12 (garnet, Ia-3d) has 120 partially occupied Li sites in its cubic cell, on `24d` and `96h`, for 56 Li atoms. Picking the 28 occupied sites out of the 60 of the primitive cell at once would mean ~10<sup>17</sup> combinations. Sequential doping makes this feasible: the Li atoms are removed one at a time, and only the best configurations of each step are continued.

```bash
pip install mace-torch
./run.sh                      # MACE-MP-0 (small), single-point energies
ENGINE="-pot LJ" ./run.sh     # only to see the workflow (~30 s; the LJ energies mean nothing)
```

## The sites

```text
$ config_enum -i LLZO.cif -prim -dry -dt Va -nd 1
Partially occupied sites in LLZO.cif (all sites are treated as occupied; use vacancies, -dt Va, to remove atoms):
  Li: 96 sites with occupancy 0.28 (26.9 atoms; 69.1 to remove)
  Li: 24 sites with occupancy 0.37 (8.9 atoms; 15.1 to remove)
Finding the primitive cell using SPGlib as requested.
48 symmetry operations were detected. 
16a: ['Zr12', ...]
24c: ['La0', ...]
24d: ['Li68', ..., 'Li79']
96h: ['Li80', ..., 'Li127', 'O20', ...]
```

With `-prim` the sites keep the labels of the cubic cell (`24d`, `96h`), although the primitive cell has half as many atoms. The refined occupancies add up to less Li than the formula (a known limitation of X-ray refinements for Li), so the target here is the nominal composition.

## What it does (`input.toml`)

- `noDopants = 32`, `seq = 1`: 32 steps, each removing one Li from the structures picked in the step before.
- `seqpick = 2`: the two configurations with the lowest reduced energy are continued, so each step has two input structures (`input_0`, `input_1`) and at most 2 × (number of Li sites left) configurations. Structures that are equivalent to one another are merged, and the same structure is not picked twice.
- The last step gives La<sub>12</sub>Li<sub>28</sub>O<sub>48</sub>Zr<sub>8</sub>; the picked structures are in `Step32/Picked/`.

The run makes about 2500 configurations in total (2, 62, 114, 113, ... per step with the LJ energies; the numbers depend on which configurations are picked, and so on the potential).

## Choices to think about

- **Single points or relaxations.** Neighbouring `24d` and `96h` sites are too close (≈1.6 Å) to be occupied together, so the first steps mainly remove such pairs. Single-point energies are used here to keep the run short; add `-opt` for relaxed energies.
- **`seqpick`**: more configurations per step explore more of the configuration space, at a cost proportional to it.
- **`seq`**: removing two Li per step (`-seq 2`) considers pairs, but gives ~1700 configurations per input structure once the symmetry is lost.
- **`-stype S`** picks the most symmetric configurations instead of the lowest in energy.
- **Final check with DFT**: run VASP on the structures in `Step32/Picked/`, or repeat the last steps with `-vasp`.

The example of the original script, up to two vacancies on the `96h` sites with one folder per configuration, still works (without the CASM files, which are no longer written):

```bash
config_enum -i LLZO.cif -prim -at Li -wl 96h -dt Va -nd 0:2 -sf -ot vasp -o POSCAR
```

It gives 1 + 1 + 33 configurations.

## Notes

- No reference energies are given: this tutorial could not be run with MACE when it was written.
- The structure is built from the data of Geiger *et al.*, Inorg. Chem. **50**, 1089 (2011).
