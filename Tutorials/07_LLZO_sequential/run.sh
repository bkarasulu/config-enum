#!/bin/bash
# Tutorial 7: LLZO by sequential doping. Equivalent to:
# config_enum -i LLZO.cif -prim -at Li -dt Va -nd 32 -seq 1 -spick 2 --mlip mace --mlip_model small
# For a quick look at the workflow only (~30 s): ENGINE="-pot LJ" ./run.sh
cd "$(dirname "$0")"
ENGINE=${ENGINE:-"--mlip mace --mlip_model small"}
config_enum --input input.toml $ENGINE -ow "$@" | tee config_enum.log
