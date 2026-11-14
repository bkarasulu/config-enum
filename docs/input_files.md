# Input files

All options can be given in a file instead of on the command line:

```bash
config_enum --input run.toml
config_enum --input run.toml -nd 3 -ow        # command-line options override the file
```

The format is chosen by the extension: TOML (`.toml`, or any other extension), JSON (`.json`) or YAML (`.yaml`, `.yml`; needs `pip install pyyaml`).

## Keys

A key is the long option name without the dashes (`noDopants`, `aType`, `seqpick`) or the short one (`nd`, `at`, `spick`); see the [Command-line reference](cli.md). Flags take `true`/`false`. Options that take several values are lists.

```toml
# Au in Cu, sequential doping
[structure]              # section names are free; they only organise the file
inpf = "Cu.cif"
rep = [2, 2, 2]

[enumeration]
aType = ["Cu"]
dType = "Au"
noDopants = 4            # or "0:2", or [1, 2]
seq = 1
seqpick = 2

[calculation]
mlip = "mace"
mlip_model = "small"
mlip_kwargs = {default_dtype = "float32"}
optimise = true
```

Unknown keys are reported as an error.

## Saving the options of a run

`--save_input run.toml` writes the options that differ from the defaults (and the required ones) to a file that can be used with `--input`.
