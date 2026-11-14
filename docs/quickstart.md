# Quick start

## 1. Look at the sites

```bash
config_enum -i LiCoO2.cif -r 2 2 1 -dry -dt Va -nd 1
```

```text
144 symmetry operations were detected. 
Space group= R-3m (166) with tol=1.00e-03
The following equivalent sites were extracted from LiCoO2.cif: 
12a: ['Co0', 'Co1', ..., 'Co11']
12b: ['Li12', 'Li13', ..., 'Li23']
24c: ['O24', 'O25', ..., 'O47']
```

The labels are the number of equivalent atoms in the cell and the Wyckoff letter. The atom indices start at 0 here; the `-nl` and `-excl` options count from 1.

## 2. Choose the sites and the dopant

| What | Option | Example |
|---|---|---|
| Sites by atom type | `-at` | `-at Li` |
| Sites by Wyckoff label | `-wl` | `-wl 12b` (with `-at`: both must match) |
| Sites by atom number (from 1) | `-nl` | `-nl 1,2,5,8:12` |
| Leave some atoms out | `-excl` | `-excl 1:4 24` |
| Dopant, or vacancy | `-dt` | `-dt Na`, `-dt Va` |
| How many | `-nd` | `-nd 2`, `-nd 0:3` |

## 3. Choose how the energies are computed

| Engine | Options |
|---|---|
| Simple ASE potential (tests only) | default: `-pot LJ`; `-pot EMT` for some metals |
| Machine-learning potential | `--mlip mace --mlip_model medium-mpa-0 -opt` |
| VASP | `-vasp -id INPUTS -np 64 -nt 2` |
| CASTEP | `-castep -iseed myseed` |

```bash
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 2 --mlip mace -opt
```

The summary is printed and written to `Step1/degeneracies`; the structures are in `Step1/`. See [Workflow and output files](workflow.md).

## 4. Many dopants: sequential doping

The number of combinations grows quickly with the number of dopants. With `-seq N`, the dopants are added N at a time, and only the best `-spick` configurations of each step are continued:

```bash
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 6 -seq 1 -spick 2 --mlip mace -opt
```

## 5. Keep the options in a file

```bash
config_enum -i Cu.cif -r 2 2 2 -at Cu -dt Au -nd 2 -pot EMT --save_input run.toml
config_enum --input run.toml -ow
```

See [Input files](input_files.md).
