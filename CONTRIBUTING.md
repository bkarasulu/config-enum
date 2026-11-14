# Contributing to config_enum

Thank you for helping to improve config_enum! Contributions of any size are welcome: bug reports, tests, examples, documentation and code.

## Working with the repository

Please always work from a clone of the GitHub repository rather than from private copies of the code. This keeps everyone's changes under version control.

```bash
git clone https://github.com/bkarasulu/config-enum.git
cd config-enum
pip install -e ".[test,lint,docs]"
git checkout -b my-feature        # new work goes on a branch
```

## Before opening a pull request

```bash
python -m pytest            # all tests must pass (~10 s; no DFT code needed)
ruff check .                # lint
python docs/gen_cli.py      # if you changed command-line options
mkdocs build --strict       # if you changed the documentation
```

- Add a test for every bug fix and every new feature (see `tests/` and `docs/development.md`).
- Keep the existing command-line options and function names working: other people's scripts use them.
- Describe the change, and how you tested it, in the pull request.

## Reporting bugs and requesting features

For bugs, please include:

- the full command line or input file;
- the input structures (or how to obtain them);
- the config_enum version (`config_enum --version`) and the code used (VASP, CASTEP, MLIP, ...);
- the error message or the unexpected output.

A small input that reproduces the problem with `-pot EMT` or `-pot LJ` is the most helpful.

## Testing new and old features (checklist for testers)

- [ ] Redo a system you studied with the original script and compare the results. The original script did not merge all symmetry-equivalent configurations (see `docs/known_issues.md`), so expect fewer configurations, with larger degeneracies, and the same energies for the configurations that remain.
- [ ] Check the sites with `-dry` and the VASP inputs of one configuration before submitting real runs.
- [ ] Run the tutorials in `Tutorials/`.
- [ ] Report anything unexpected as an issue.

## Licence

config_enum is released under the GNU General Public License v3.0 or later (see `LICENSE`). By contributing, you agree that your contributions are released under the same licence.
