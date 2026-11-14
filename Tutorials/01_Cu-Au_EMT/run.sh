#!/bin/bash
# Tutorial 1: Au in Cu with EMT. Equivalent to:
# config_enum -i Cu.cif -r 2 2 2 -at Cu -dt Au -nd 4 -seq 1 -spick 2 -stype E -pot EMT -opt --fmax 0.02 -t 300
cd "$(dirname "$0")"
config_enum --input input.toml -ow "$@" | tee config_enum.log
