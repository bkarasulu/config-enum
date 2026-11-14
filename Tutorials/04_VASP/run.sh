#!/bin/bash
# Tutorial 4: Li vacancies in LiCoO2 with VASP (sequential doping). Run it in a job script on your cluster.
# Needs VASP_PP_PATH (for the POTCARs) and either VASP_COMMAND or the -np/-nt/-mpi options.
cd "$(dirname "$0")"
config_enum -i LiCoO2.cif -r 2 2 1 -at Li -dt Va -nd 4 -seq 1 -spick 2 \
            -vasp -id INPUTS -xc pbe -hubU Co d 3.32 -mgm Co 1.0 -np 64 -nt 2 -mpi -rs "$@" | tee config_enum.log
