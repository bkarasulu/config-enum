# Tutorial 4: Li vacancies in LiCoO2 with VASP

This tutorial shows how to set up a sequential enumeration with VASP. It is meant to be run on a cluster; the settings in `INPUTS/` are examples and must be converged for your system.

```bash
./run.sh            # inside a job script
```

## Set-up

- `INPUTS/INCAR` and `INPUTS/KPOINTS` are copied to the run folder of every configuration (`-id INPUTS`). The INCAR decides what is done (here a relaxation of the atoms, `IBRION = 2`, `ISIF = 2`).
- The POTCARs are made by ASE from `$VASP_PP_PATH` (the recommended potentials for `-xc pbe`). Put a `POTCAR` in `INPUTS/` only for reference; ASE writes its own.
- `-hubU Co d 3.32` and `-mgm Co 1.0` add the Hubbard U (LDAU tags) and the initial magnetic moments to the INCAR.
- `-np 64 -nt 2 -mpi` runs two configurations at a time, each with `mpirun -np 32 vasp_std`. Without `-mpi`, `srun` is used (with `-nn` nodes). `VASP_COMMAND`, if set, is used instead.

## What it does

`-nd 4 -seq 1 -spick 2`: four Li vacancies are made one at a time. At each step, every symmetrically-unique configuration is relaxed with VASP in `StepN/input_M/config_K/`, and the two configurations with the lowest reduced energy are the starting structures of the next step.

## Restarts

With `-rs`, an interrupted run can be started again with the same command: finished VASP runs are read in, and unfinished ones are continued from their `CONTCAR` (the files of the previous attempt are kept in `RUN1/`, `RUN2/`, ...). Use `-ow` to start from scratch (this deletes the `Step*` folders).

A configuration whose calculation fails is skipped with a warning and gets an energy of 0 as a placeholder (it is then never picked while other configurations are available); with `-noskip` the run stops instead.

## Check before you submit

```bash
config_enum -i LiCoO2.cif -r 2 2 1 -dry -dt Va -nd 1          # the sites
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 4 -seq 1 -spick 2 -pot LJ   # number of configurations per step (seconds)
```
