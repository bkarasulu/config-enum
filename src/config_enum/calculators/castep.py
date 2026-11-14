"""Quick CASTEP calculations through ASE."""
import glob
import os
import shutil

from ..log import say, warn


def castep_command(nprocs=32, ntasks=1, nnodes=1, mpirun=False, castepexe='castep.mpi'):
    """The command to run CASTEP: $CASTEP_COMMAND, or srun/mpirun with the cores shared by the tasks."""
    cascmd = os.getenv("CASTEP_COMMAND")
    if not cascmd:
        if mpirun:
            cascmd = "mpirun -np %d %s" % (nprocs / ntasks, castepexe)
        else:
            cascmd = 'srun -n %d -N %d %s' % (nprocs / ntasks, max(1, nnodes / ntasks), castepexe)
    return cascmd


def call_castep(atoms, label='config', directory='CASTEP-tmp', cascmd=None, inpseed=None, base_dir='.', pp_path=None):
    """Run a quick, low-accuracy CASTEP single-point calculation (or with the settings of <inpseed>.param/.cell).

    Args:
        atoms: structure.
        label: seed name of the calculation.
        directory: run folder.
        cascmd: CASTEP command.
        inpseed: seed of the .param and .cell files (in base_dir) with the settings to use (keywords, pseudopotentials
            and k-points).
        base_dir: folder with the <inpseed> files and the .usp pseudopotential files.
        pp_path: folder of the pseudopotentials (default: $CASTEP_PP_PATH or $HOME/APPS/pspots).

    Returns (energy, atoms, pressure in GPa)
    """
    import ase.calculators.castep
    import ase.io.castep

    atoms.calc = None
    calc = ase.calculators.castep.Castep()

    if cascmd:
        calc._castep_command = cascmd
    calc._castep_pp_path = pp_path or os.environ.get('CASTEP_PP_PATH') or os.path.join(os.path.expanduser('~'), "APPS", "pspots")
    calc._directory = directory
    calc.cell.kpoint_mp_grid = '1 1 1'
    calc.param.xc_functional = 'PBE'
    calc.param.cut_off_energy = 300
    calc.param.max_scf_cycles = 100
    calc.param.mix_charge_amp = 0.1  # for treating slow convergence
    calc.param.task = 'SinglePoint'
    calc.param.write_checkpoint = 'none'
    calc._rename_existing_dir = False
    calc._copy_pspots = True
    calc._label = label
    atoms.calc = calc

    if inpseed:
        # read the block species_pot and kpoint_mp_grid info from the ref cell file!!
        seed = os.path.join(base_dir, inpseed)
        say('Reading CASTEP params from %s' % seed)
        atoms.calc.merge_param('%s.param' % seed)
        atoms2 = ase.io.castep.read_castep_cell('%s.cell' % seed)

        try:
            val = atoms2.calc.cell.species_pot.value
            val = [x.split() for x in val.split('\n')]
            val.pop(0)
            atoms.calc.cell.species_pot = val
        except Exception:
            warn('Cannot assign the pspots defined in the reference cell file.')
        try:
            atoms.calc.cell.kpoints_mp_grid = atoms2.calc.cell.kpoints_mp_grid.value
        except Exception:
            pass
        try:
            atoms.calc.cell.kpoints_mp_spacing = atoms2.calc.cell.kpoints_mp_spacing.value
        except Exception:
            pass
        try:
            atoms.calc.cell.symmetry_generate = True
        except Exception:
            pass
        try:
            atoms.calc.cell.symmetry_tol = atoms2.calc.cell.symmetry_tol.value
        except Exception:
            atoms.calc.cell.symmetry_tol = 0.00001
        atoms.calc.cell.snap_to_symmetry = True

    os.makedirs(directory, exist_ok=True)
    for fn in glob.glob(os.path.join(base_dir, '*.usp')):
        shutil.copy(fn, directory)
    e2 = atoms.get_potential_energy()
    if e2 is None:
        e2 = 0.0
    try:
        P = atoms.info['pressure']
    except KeyError:
        P = 0.
    return e2, atoms, P
