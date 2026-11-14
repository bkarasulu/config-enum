"""VASP calculations through ASE (inputs are taken from the INPUTS folder, see -id)."""
import glob
import os
import re
import shutil
from os import getenv

import ase.io

from ..log import say, warn


def handle_magmoms(atoms, magmoms):
    """
    The magmoms variable should be a list parsed by argparse in the form:
    [...] -mgm Fe 5 Nb 0.6 O 0.6 [...]
    which is then converted to a dictionary:
    d = {
        'Fe': 5.,
        'Nb': 0.6,
        'O': 0.6
        }
    """
    if magmoms is None or len(magmoms) == 0:
        return atoms
    else:
        say("Magnetic moments are set with ASE: %s" % magmoms)

    elements = magmoms[::2]
    values = [float(v) for v in magmoms[1::2]]
    d = dict(zip(elements, values))
    init_mgm = []
    for s in atoms.symbols:
        if s not in elements:
            init_mgm.append(0)
        else:
            init_mgm.append(d[s])
    atoms.set_initial_magnetic_moments(init_mgm)
    return atoms


def set_hubU(atoms, hubU):
    """Set the Hubbard U parameters of the VASP calculator of atoms; hubU = {'Fe': {'d': 2.7}, 'Sm': {'f': 6.1}}"""
    if len(hubU) > 0:
        atoms.calc.set(ldau=True, ldauprint=0, ldautype=2, lmaxmix=6)
        ats = sorted(hubU.keys())
        ldauu = []
        ldaul = []
        asyms = list(dict.fromkeys(atoms.get_chemical_symbols()))  # atom types in the order of appearance
        for at in asyms:
            if at in ats:
                if len(hubU[at].keys()) > 1:
                    warn("WARNING: cannot set multiple U values for the same atom type (%s), setting to 0... " % at)
                    ldauu.append(0)
                    ldaul.append(0)
                    continue
                orb = "".join(hubU[at].keys())
                ldauu.append(float(hubU[at][orb]))
                if orb == "p":
                    ldaul.append(1)
                elif orb == "d":
                    ldaul.append(2)
                elif orb == "f":
                    ldaul.append(3)
                else:
                    raise ValueError("Hubbard U: unknown orbital type '%s' for %s (use p, d or f)" % (orb, at))
            else:
                ldaul.append(0)
                ldauu.append(0)
        atoms.calc.set(ldauj=[0 for i in range(len(asyms))])
        atoms.calc.set(ldauu=ldauu, ldaul=ldaul)
        # This is essential for accurate total energies and band structure calculations for f-elements (e.g. ceria),
        # all 3d-elements (transition metal oxides), and magnetic atoms in the 2nd row (B-F atom), in particular if
        # LDA+U or hybrid functionals or meta-GGAs are used, since these functionals often result in aspherical
        # charge densities.
        atoms.calc.set(lasph=True)
    return atoms


def vasp_continue(ifDel=0, wdir="."):
    """Store the data of the previous run in a RUN<n> folder and continue from CONTCAR."""
    ind = 0
    for fn in os.listdir(wdir):
        if fn[0:3] == "RUN" and fn[3:].isdigit():
            ind = max(ind, int(fn[3:]))

    run = "RUN%d" % (ind + 1)
    say("Previous data stored in %s." % run)
    os.makedirs(os.path.join(wdir, run))
    for fn in os.listdir(wdir):
        if os.path.isfile(os.path.join(wdir, fn)):
            shutil.copy(os.path.join(wdir, fn), os.path.join(wdir, run, fn))
    contcar = os.path.join(wdir, "CONTCAR")
    if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
        shutil.copy(contcar, os.path.join(wdir, "POSCAR"))
    if ifDel:
        for fn in ("WAVECAR", "CHGCAR"):
            if os.path.exists(os.path.join(wdir, fn)):
                os.remove(os.path.join(wdir, fn))


def get_potcar(elements, xc, vasppp=None):
    """The POTCAR text for the given elements from $VASP_PP_PATH/<xc>/<element>/POTCAR."""
    if vasppp is None:
        vasppp = os.environ['VASP_PP_PATH']
    ppdir = os.path.join(vasppp, xc)

    pp = ''
    for e in elements:
        with open(os.path.join(ppdir, e, 'POTCAR')) as f:
            pp += f.read()

    return pp


def make_potcar(xc="potpaw_PBE", elements=None, wDir='./'):
    """Compile the POTCAR in wDir for the elements (read from the POSCAR in wDir if not given)."""
    try:
        vasppp = os.environ['VASP_PP_PATH']
    except KeyError:
        raise RuntimeError("VASP_PP_PATH is not defined, but needed for VASP calculation setup.")

    say('Making potcar using %s' % vasppp)

    if elements is None:  # Try to read it from POSCAR if not explicitly given.
        try:  # VASP5 format.
            with open(os.path.join(wDir, "POSCAR")) as f:
                elements = f.readlines()[5].strip().split()
        except (OSError, IndexError):
            raise RuntimeError('Elements not given and POSCAR not found in %s' % wDir)

    pp = get_potcar(elements, xc, vasppp=vasppp)

    with open(os.path.join(wDir, "POTCAR"), 'w') as p:
        p.write(pp)


def vasp_command(exe=None, nprocs=32, ntasks=1, nnodes=1, mpirun=False, vaspexe='vasp_std'):
    """The command to run VASP: the argument, $VASP_COMMAND, or srun/mpirun with the cores shared by the tasks."""
    if not exe:
        exe = getenv('VASP_COMMAND')
    if not exe:
        if mpirun:
            exe = 'mpirun -np %d %s' % (nprocs / ntasks, vaspexe)
        else:
            exe = 'srun -n %d -N %d %s' % (nprocs / ntasks, max(1, nnodes / ntasks), vaspexe)
    return exe


def read_pressure(outcar):
    """External pressure + Pulay stress of the last step in an OUTCAR (GPa); ASE does not get it right when restarting."""
    P = 0.
    try:
        with open(outcar, errors="replace") as f:
            for ln in f:
                if 'pressure' in ln:
                    m = re.search(r'external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)', ln)
                    if m:
                        P = (float(m.group(1)) + float(m.group(2))) * 0.1  # kB to GPa
    except (OSError, ValueError):
        pass
    return P


def copy_vasp_inputs(idir, wdir):
    """Copy INCAR, KPOINTS and POTCAR (if they exist) from the input folder to the run folder."""
    for fn in ('INCAR', 'POTCAR', 'KPOINTS'):
        src = os.path.join(idir, fn)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(wdir, fn))


def prepare_vasp(wdir, xc='pbe', mgms=None, hubU=None):
    """Set up a new (or continuation) VASP run in wdir from the INCAR, KPOINTS and POSCAR there.

    Returns the atoms with the ASE VASP calculator attached.
    """
    from ase.calculators.vasp import Vasp

    hubU = hubU or {}
    calc = Vasp(directory=wdir)
    incar = os.path.join(wdir, 'INCAR')
    if os.path.exists(incar):
        calc.read_incar(filename=incar)
    else:
        warn("WARNING: no INCAR in %s; the default settings of ASE are used." % wdir)
    kpoints = os.path.join(wdir, 'KPOINTS')
    if os.path.exists(kpoints):  # otherwise ASE would overwrite it with a Gamma-point only KPOINTS
        try:
            calc.read_kpoints(filename=kpoints)
        except Exception as err:
            warn("WARNING: could not read %s (%s)" % (kpoints, err))
    if os.path.exists(os.path.join(wdir, "OUTCAR")):
        vasp_continue(wdir=wdir)
    atoms = ase.io.read(os.path.join(wdir, "POSCAR"), format="vasp")
    calc.set(xc=xc)
    calc.set(setups='recommended')
    atoms.calc = calc

    # Adding MAGMOM and ISPIN to INCAR if -mgm or --magmoms is defined.
    atoms = handle_magmoms(atoms=atoms, magmoms=mgms)
    if len(hubU) > 0:
        atoms = set_hubU(atoms, hubU)
    return atoms


def call_vasp_v2(wdir='.', exe=None, xc='pbe', mgms=None, hubU=None, makepotcar=False, potcarxc="potpaw_PBE"):
    """Run VASP in the folder wdir (which contains POSCAR, INCAR and KPOINTS), or read the results if it is already
    converged; a run that is not converged is continued. Returns (energy, atoms)."""
    from ase.calculators.vasp import Vasp

    if not exe:
        exe = getenv('VASP_COMMAND')
    if not exe:
        exe = 'vasp_std'

    os.environ['VASP_COMMAND'] = exe
    say("VASP command: ", exe)
    say('Working dir: %s' % wdir)

    if makepotcar:
        make_potcar(xc=potcarxc, wDir=wdir)

    flag = 0  # whether to start a new/continuation run
    atoms = None
    try:
        calc = Vasp(restart=True, directory=wdir)
        atoms = calc.get_atoms()
        say("\nVASP run was read successfully from OUTCAR.")
        if calc.read_convergence():
            say('Geom opt already converged...')
            flag = 0
        else:
            say('Geom opt not converged; running a continuation job...')
            flag = 1
    except Exception:
        say("VASP run could not be read, starting a new run...")
        flag = 1

    if flag:
        atoms = prepare_vasp(wdir, xc=xc, mgms=mgms, hubU=hubU)

    Ep = atoms.get_potential_energy()
    return Ep, atoms


def clean_vasp_files(wdir, patterns=("WAVECAR", "CHGCAR", "CHG")):
    """Delete large output files in a run folder."""
    for pat in patterns:
        for fn in glob.glob(os.path.join(wdir, pat)):
            os.remove(fn)
