#!/bin/bash
# Tutorial 5: point defects in Li2S (antifluorite). Three runs, each in its own folder.
# The energies are computed with MACE by default; for a quick look at the configurations only: ENGINE="-pot LJ" ./run.sh
cd "$(dirname "$0")"
ENGINE=${ENGINE:-"--mlip mace --mlip_model small -opt"}

# (a) two Li vacancies in the 2x2x2 supercell (96 atoms)
mkdir -p a_vacancies && (cd a_vacancies && config_enum -i ../Li2S.cif -r 2 2 2 -at Li -dt Va -nd 2 $ENGINE -ow "$@" | tee config_enum.log)

# (b) Cl on a S site, with the nearest Li removed for charge balance (Cl_S + V_Li)
mkdir -p b_Cl_doping && (cd b_Cl_doping && config_enum -i ../Li2S.cif -r 2 2 2 -at S -dt Cl -nd 1 \
    -mod "remove Li nearest 1 of Cl" $ENGINE -ow "$@" | tee config_enum.log)

# (c) Li Frenkel pair: Li on an octahedral interstice (4b; 32b in the supercell), and the nearest lattice Li removed
mkdir -p c_frenkel && (cd c_frenkel && config_enum -i ../Li2S.cif -ai -r 2 2 2 -at X -wl 32b -dt Li -nd 1 \
    -mod "remove Li nearest 1 of dopant" $ENGINE -ow "$@" | tee config_enum.log)
