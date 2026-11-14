# DFT codes

## VASP (`-vasp`)

Each configuration is run in its own folder, `StepN/input_M/config_K/`:

1. `POSCAR` is written, and `INCAR` and `KPOINTS` (and `POTCAR`, if present) are copied from the input folder (`-id`, default `INPUTS/`, relative to the folder where config_enum is started).
2. ASE's VASP calculator is set up from the `INCAR` and `KPOINTS`. It writes the `POTCAR` from `$VASP_PP_PATH` with the recommended potentials for `-xc` (default `pbe`), and adds
    - the initial magnetic moments of `-mgm`, e.g. `-mgm Fe 5.0 O 0.6` (atom types that are not given get 0);
    - the Hubbard U of `-hubU`, e.g. `-hubU Fe d 2.7, Sm f 6.1` (LDAU, LDAUTYPE = 2, LDAUL/LDAUU in the order of the POSCAR, LASPH).
3. VASP is run. What is computed (single point, relaxation, cell relaxation) is decided by your `INCAR`.

The command is `$VASP_COMMAND` if set; otherwise it is built from the options: `srun -n <np/nt> -N <nn/nt> <vexe>`, or `mpirun -np <np/nt> <vexe>` with `-mpi`. `-nt` configurations are run at the same time, sharing the `-np` cores.

`-mp` compiles the POTCAR from `$VASP_PP_PATH/<--potcarxc>/<element>/POTCAR` in the order of the POSCAR, for the folders written with `-sf -copy` and before a VASP run.

### Restarts

With `-rs`, the same command can be run again: a finished (converged) run is read in and not redone, and an unconverged run is continued from its `CONTCAR`. The files of the earlier attempt are kept in `RUN1/`, `RUN2/`, ... The pressure is read from the `OUTCAR`.

### Failed calculations

A configuration whose calculation fails is skipped with a warning; it is listed with an energy of 0 and space group P1, and it is not picked while other configurations are available. With `-noskip` the run stops at the first failure.

## CASTEP (`-castep`)

Without further options, a quick, low-accuracy single-point calculation is done (PBE, 300 eV cut-off, Γ point) to estimate the energies. With `-iseed <seed>`, the settings of `<seed>.param` and the pseudopotentials, k-points and symmetry tolerance of `<seed>.cell` are used instead (both in the folder where config_enum is started). `.usp` files in that folder are copied to the run folder (`StepN/input_M/CASTEP-tmp/`).

The command is `$CASTEP_COMMAND`, or built from `-np/-nt/-nn/-mpi/-cexe` like the VASP one. The pseudopotential folder is `$CASTEP_PP_PATH` (default `$HOME/APPS/pspots`, as in the original script).

!!! note
    ASE's CASTEP interface needs a working CASTEP installation to set itself up (it creates `castep_keywords.json` on first use).

## Writing inputs only

To prepare input folders without running anything, use the default (LJ) energies and write the structures in separate folders:

```bash
config_enum -i LiCoO2.cif -at Li -dt Va -nd 2 -sf -ot vasp -o POSCAR -copy -id INPUTS -mp
```

This gives `Step1/input_0/<K>/POSCAR` with the files of `INPUTS/` (and a POTCAR with `-mp`) for each configuration.
