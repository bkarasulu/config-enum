# Tutorial 1: Au in Cu with EMT (seconds)

This tutorial dopes a 32-atom fcc Cu cell with four Au atoms, one at a time (sequential doping), using ASE's EMT potential. It needs no DFT code or ML package and takes about 2 s.

```bash
./run.sh            # = config_enum --input input.toml -ow
```

## What it does (`input.toml`)

1. Builds the 2×2×2 supercell of the cubic cell (`rep`). All 32 Cu sites are equivalent (Wyckoff site `32a` of the supercell) and are selected with `aType = ["Cu"]`.
2. **Step 1**: one Au atom. There are 32 ways, all equivalent: one configuration with degeneracy 32.
3. **Step 2**: a second Au atom is added to the configuration picked in Step 1. The 31 remaining sites give 5 symmetrically-unique configurations (first to fifth neighbours in this cell). Each is relaxed with EMT (`optimise`), and the two with the lowest reduced energy are picked (`seqpick = 2`).
4. **Steps 3 and 4**: the same for each of the picked configurations (`input_0`, `input_1`). Configurations that are equivalent to another one from the same input structure are merged, and equivalent configurations made from different input structures are not picked twice.

## Expected results

`Step2/degeneracies` (also printed on the screen):

```text
Step1/input_0/config_0.res
Config-ID   Degens      E_EMT       E_reduced  Prob.    Space Group
--------------------------------------------------------------------------------
 config_0     12        -0.10295      -0.16719  0.043    Cmmm (65) 
 config_1     3         -0.19826      -0.22666  0.430    P4/mmm (123) 
 config_2     12        -0.16373      -0.22796  0.452    Cmcm (63) 
 config_3     3         -0.13619      -0.16459  0.039    P4/mmm (123) 
 config_4     1         -0.16204      -0.16204  0.035    Im-3m (229) 
...
Picked configs for Step 2: 
Step2/input_0/config_2.res  Step2/input_0/config_1.res
```

- `Degens` is the number of equivalent ways of making the configuration from the input structure (12 + 3 + 12 + 3 + 1 = 31 sites).
- `E_reduced` = E − k<sub>B</sub>T ln(degeneracy) and `Prob.` are the reduced energy and the probability at `temp` (see `docs/theory.md`).
- Nearest-neighbour Au pairs (`config_0`) are the least favourable: with EMT, Au atoms in Cu repel each other at short distance.

The numbers of configurations are 1, 5, 13 and 24 for the four steps.

## Things to try

- `./run.sh -stype S`: pick the configurations with the highest symmetry instead of the lowest energy.
- `./run.sh -seq 2`: add two Au atoms per step (more configurations per step, fewer steps).
- `./run.sh -seq 0 -nd 2`: no sequential doping, all pairs at once (the same five configurations as in Step 2, with the degeneracies of the pairs: 192, 48, 192, 48, 16).
- `./run.sh -dt Va`: vacancies instead of Au.
- Look at the structures: `ase gui Step4/Picked/*.res`; the relaxed ones are `Step4/Au*-config_enum-*.res`.
