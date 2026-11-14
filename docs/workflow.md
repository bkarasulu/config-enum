# Workflow and output files

## Steps

Without `-seq` there is one step (`Step1/`), in which all the dopants of `-nd` are placed at once; several numbers can be given (`-nd 0:3`). With `-seq N`, the dopants are added N at a time until `-nd` is reached (a single number), in steps `Step1/`, `Step2/`, ...

For every input structure of a step (the files of `-i` in Step 1, the picked configurations afterwards):

1. **Structure** (Step 1 only): the interstices are added (`-ai`), the primitive cell is found (`-prim`) and the supercell is built (`-r`), in this order.
2. **Symmetry**: the symmetry operations and the equivalent sites are found with spglib (`-tol`). In Step 1 they are those of the structure itself, or of the reference structure (`-ref`). In later steps they are those of the doped structure (see [Theory](theory.md#sequential-doping) and `-ssym`). `-nosymm` switches symmetry off; `-dry` stops here.
3. **Selection** (Step 1 only): `-at`, `-wl` and `-nl` select the sites, `-excl` removes atoms from the selection. In later steps the sites that are left are used.
4. **Enumeration**: all combinations of the dopants on the selected sites; one configuration of each set of equivalent ones is kept, with its degeneracy. If several atom types are selected (`-at Li Nb`), only configurations in which all of them are doped are kept.
5. **Doping**: the selected atoms are replaced by `-dt`, or deleted for `Va`. With `-dexcl`, the atoms of `-excl` are doped too, in every configuration, with `-dexcl_type` (default: `-dt`). Then the modifications (`-mod`) are applied and the unused interstices are removed.
6. **Final check** (`-ase`): the whole structures are compared with ASE's `SymmetryEquivalenceCheck`, and equivalent ones are merged.
7. **Output** of the structures (`-ot`, `-o`, `-sf`, `-copy`).
8. **Energies** with the chosen engine, `-nt` configurations at a time.
9. **Thermodynamics**: reduced energies and probabilities at `-t` ([Theory](theory.md#configurational-thermodynamics)).

At the end of the step the configurations of all its input structures are sorted (`-stype`) and the first `-spick` are picked.

## Modifications (`-mod`)

```text
-mod '<replace|remove> <atom type|all> <within|nearest> <value> of <reference> [with <new type>]'
```

- `within 2.5`: all atoms of the given type within 2.5 Å of a reference atom; `nearest 2`: the two nearest ones of each reference atom.
- The reference is an atom type, or `dopant` for the doped (or vacant) sites of the configuration.
- `replace ... with Va` and `remove` delete the atoms.
- Several statements are separated by `;` and applied in order.
- If several atoms are at the same distance, `nearest` takes them in the order of the structure; when they are not equivalent by symmetry, the choice is arbitrary. In that case enumerate the second defect in a separate run instead.

## Partially occupied sites

Structures with partial occupancies (e.g. CIF files from a database) are read with all sites occupied. config_enum lists these sites and the number of atoms to remove; make the ordered models with `-dt Va` on those sites (see Tutorials 6 and 7).

Examples: `-mod 'remove Li nearest 1 of Al'` (one Li removed for each Al), `-mod 'replace O within 2.2 of dopant with F'`.

## Interstitial sites

With `-ai -at X`, the unoccupied Wyckoff positions of the space group are added as `X` atoms before the enumeration, and `-dt` is placed on them. They appear in the list of sites and can be selected with `-wl`. The `X` atoms that are not used are removed from the output; in a sequential doping they stay available for the next steps.

```bash
config_enum -i Cu.cif -ai -at X -dt H -nd 1 -pot EMT        # H on the octahedral (4b), tetrahedral (8c) and 24d sites of fcc Cu
config_enum -i Cu.cif -ai -at X -wl 4b -dt H -nd 2 -seq 1   # octahedral sites only, two H atoms one after the other
```

Only the positions without free parameters are added (e.g. `4b`, `8c` and `24d` in Fm-3m, `8b`, `16c` and `16d` in diamond). A position such as (x,x,x) is a line of possible sites; add atoms on such sites to the input structure yourself (as `X` atoms, selected with `-at X` without `-ai`).

The positions of all space groups and settings are part of the code (`config_enum/wyckoff_positions.py`). They are those of *International Tables for Crystallography, Vol. A, 5th ed., edited by Th. Hahn (Kluwer, Dordrecht, 2005)*, taken from the `wyckoff.dat` file of VESTA that the original script read:

```python
from config_enum.symmetry import wyckoff_positions
from config_enum.wyckoff_positions import WYCKOFF_POSITIONS

WYCKOFF_POSITIONS[(225, 1)]             # ('F m -3 m', (('192l', '1', '(x,y,z) ...'), ..., ('8c', '-43m', '(1/4,1/4,1/4) (1/4,1/4,3/4)'), ...))
wyckoff_positions(225, fixed_only=True) # {'225-1': {'24d': [[0.0, 0.25, 0.25], ...], '8c': [...], '4b': [...], '4a': [...]}}
```

Other data can be used with `--wyckoff_file FILE` (or `$CONFIG_ENUM_WYCKOFF`), in the format of the original file: for each space-group setting a header line starting with the space-group number and the setting, then for each Wyckoff position a line `<multiplicity> <letter> <site symmetry>` followed by lines of coordinates in brackets; settings are separated by empty lines.

Limitations: the structure must be a conventional cell in the first setting of its space group in the table (e.g. origin choice 1, unique axis b). The interstices are added before `-prim` and `-r` are applied.

## Output files

```text
Step1/
  degeneracies                      the summary of the step (also printed)
  input_0/                          one folder per input structure of the step
    config_0.res, config_1.res ...  the configurations: dopants on the lattice sites, not optimised
    config_0/                       VASP run folder of config_0 (with -vasp); CASTEP-tmp/ with -castep
  <formula>-config_enum-<id>.res    structures after the calculation; energy, pressure and space group in the title line
  Picked/input_0-config_3.res       copies of the picked configurations
  Structures-are-optimised          marker: the <formula>-...res files are relaxed (VASP, or -opt)
Step2/ ...
results.json
```

- `-ot` sets the format of the configurations (any ASE format, e.g. `vasp`, `cif`, `xyz`; default `res`); `-o NAME` the name (`NAME_<ID>.<ext>`).
- `-sf` writes each configuration to its own folder, `input_0/<ID>/`; with `-o POSCAR` the file is called `POSCAR`, and `-copy` copies the files of the input folder (`-id`) there.
- The structures used as the input of the next step are always the unoptimised ones.

### The summary

```text
Step1/input_0/config_0.res
Config-ID   Degens      E_EMT       E_reduced  Prob.    Space Group
--------------------------------------------------------------------------------
 config_0     12        -0.10295      -0.16719  0.043    Cmmm (65) 
 ...
All sorted configs for Step 2
        Config-ID               Energy     SG 
------------------------------------------------------------
Step2/input_0/config_2.res -0.22796... 63
 ...
Picked configs for Step 2: 
Step2/input_0/config_2.res  Step2/input_0/config_1.res
```

The energy column is `E_DFT` for VASP and CASTEP, and `E_<potential>` otherwise. The space group is that of the structure after the calculation. In the sorted list the energy is the reduced one.

### results.json

For each step: the numbers of dopants, the input structures, the picked files, and for each configuration its name, file, doped sites, degeneracy, energy, reduced energy, probability, pressure, space group, formula, result file, and whether it failed or was picked.

## Existing output

If `Step1/` exists, config_enum stops. `-ow` deletes all `Step<N>` folders and starts again; `-rs` keeps them and reuses finished VASP runs (the enumeration is repeated, which is fast, and gives the same configurations).
