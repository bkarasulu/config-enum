#!/bin/bash
# Tutorial 2: Li vacancies in LiCoO2 with MACE. Equivalent to:
# config_enum -i LiCoO2.cif -r 2 2 1 -at Li -wl 12b -dt Va -nd 1:3 --mlip mace --mlip_model small -opt
cd "$(dirname "$0")"
config_enum --input input.toml -ow "$@" | tee config_enum.log
