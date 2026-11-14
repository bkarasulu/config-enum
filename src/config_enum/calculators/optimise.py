"""Shared ASE geometry optimiser used by the ASE-based engines (classical and machine-learning potentials)."""
import warnings

from ase.optimize import BFGS, FIRE, LBFGS

from ..log import debug, verbose, warn

OPTIMISERS = {"BFGS": BFGS, "FIRE": FIRE, "LBFGS": LBFGS}


def get_optimiser(name):
    """Return the ASE optimiser class for a name (BFGS, FIRE or LBFGS, case-insensitive)."""
    try:
        return OPTIMISERS[name.upper()]
    except KeyError:
        raise ValueError("Unknown optimiser '%s'; choose from %s" % (name, ", ".join(OPTIMISERS)))


def relax(atoms, calc, fmax=0.05, steps=200, optimiser="BFGS", relax_cell=False, trajectory=None, logfile=None):
    """Relax `atoms` in place with an ASE optimiser and return (energy, atoms).

    Args:
        atoms (ase.Atoms): structure to relax (modified in place).
        calc: ASE calculator to attach.
        fmax (float): force convergence criterion (eV/A).
        steps (int): maximum number of optimisation steps.
        optimiser (str): BFGS, FIRE or LBFGS.
        relax_cell (bool): also relax the cell (FrechetCellFilter).
        trajectory (str): optional trajectory file name.
        logfile: optimiser log ('-' = stdout); default: stdout with -v.
    """
    atoms.calc = calc
    target = atoms
    if relax_cell:
        from ase.filters import FrechetCellFilter
        target = FrechetCellFilter(atoms)
    if logfile is None and verbose():
        logfile = "-"
    opt = get_optimiser(optimiser)(target, trajectory=trajectory, logfile=logfile)
    try:
        with warnings.catch_warnings():   # harmless round-off notes of scipy's logm in FrechetCellFilter
            warnings.filterwarnings("ignore", message="logm result may be inaccurate")
            converged = opt.run(fmax=fmax, steps=steps)
    except Exception as err:
        warn("relax: WARNING: optimisation failed (returning the energy of the last geometry): %s" % err)
        converged = False
    if not converged:
        warn("relax: WARNING: not converged to fmax=%.3f eV/A within %d steps" % (fmax, steps))
    E = atoms.get_potential_energy()
    debug("  %s after %d steps: E = %.5f eV" % ("converged" if converged else "stopped", opt.nsteps, E))
    return E, atoms
