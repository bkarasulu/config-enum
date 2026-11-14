#!/bin/bash
# Tutorial 6: site disorder in Li3PS4. Equivalent to:
# config_enum -i Li3PS4.cif -at Li -wl 4b 4c -dt Va -nd 4 --mlip mace --mlip_model small -opt
# For a quick look at the configurations only: ENGINE="-pot LJ" ./run.sh
cd "$(dirname "$0")"
ENGINE=${ENGINE:-"--mlip mace --mlip_model small -opt"}
config_enum --input input.toml $ENGINE -ow "$@" | tee config_enum.log
