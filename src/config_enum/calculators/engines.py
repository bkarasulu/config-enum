"""Engines: a common interface to run VASP, CASTEP and ASE calculators on a configuration."""
import os

from ase.units import GPa

from ..log import say
from .base import CalcSettings, Engine
from .castep import call_castep, castep_command
from .classical import classical_calculator
from .mlip import lazy_calculator_factory
from .optimise import relax
from .vasp import call_vasp_v2, copy_vasp_inputs, read_pressure, vasp_command


class AseEngine(Engine):
    """Any ASE calculator (classical or machine-learning potential): single point, or geometry optimisation (-opt).

    Args:
        calc_factory: a function (or class) that returns the ASE calculator, or a calculator.
        label: name used in the output.
    """

    def __init__(self, calc_factory, label="ASE", settings=None):
        super().__init__(settings)
        is_instance = not isinstance(calc_factory, type) and hasattr(calc_factory, "get_potential_energy")
        self._factory = (lambda: calc_factory) if is_instance else calc_factory
        self._calc = None
        self.label = label

    def get_calculator(self):
        if self._calc is None:
            self._calc = self._factory()
        return self._calc

    def run(self, atoms, name="config", directory="."):
        s = self.settings
        calc = self.get_calculator()
        atoms.calc = calc
        if s.optimise:
            E, atoms = relax(atoms, calc, fmax=s.fmax, steps=s.steps, optimiser=s.optimiser, relax_cell=s.relax_cell)
        else:
            E = atoms.get_potential_energy()
        P = 0.
        try:
            stress = atoms.get_stress(voigt=True)
            P = float(-(stress[0] + stress[1] + stress[2]) / 3. / GPa)
        except Exception:
            pass
        return float(E), atoms, P

    def describe(self):
        s = self.settings
        return "%s (%s)" % (self.label, "geometry optimisation, fmax=%.3f eV/A%s" % (s.fmax, ", cell relaxed" if s.relax_cell else "")
                            if s.optimise else "single point")


class VaspEngine(Engine):
    """VASP through ASE; INCAR, KPOINTS (and POTCAR) are taken from the settings.idir folder. Each configuration is run
    in <directory>/<name>; converged runs are read in and unconverged ones are continued."""
    label = "VASP"
    is_dft = True

    def __init__(self, settings=None, vaspexe="vasp_std"):
        super().__init__(settings)
        s = self.settings
        self.exe = vasp_command(s.exe, nprocs=s.nprocs, ntasks=s.ntasks, nnodes=s.nnodes, mpirun=s.mpirun, vaspexe=vaspexe)

    def run(self, atoms, name="config", directory="."):
        s = self.settings
        wdir = os.path.join(directory, name)
        os.makedirs(wdir, exist_ok=True)
        say('\nWorking on ', name)
        if not os.path.exists(os.path.join(wdir, 'OUTCAR')):  # a new run
            atoms.write(os.path.join(wdir, 'POSCAR'), format='vasp', vasp5=True)
        idir = s.idir if os.path.isabs(s.idir) else os.path.join(s.base_dir, s.idir)
        copy_vasp_inputs(idir, wdir)
        E, atoms = call_vasp_v2(wdir, exe=self.exe, xc=s.xc, mgms=s.magmoms, hubU=s.hubU, makepotcar=s.makepotcar,
                                potcarxc=s.potcarxc)
        P = read_pressure(os.path.join(wdir, 'OUTCAR'))
        return float(E), atoms, P

    def describe(self):
        return "VASP (command: %s)" % self.exe


class CastepEngine(Engine):
    """Quick CASTEP single points through ASE (or the settings of the -iseed .param/.cell files)."""
    label = "CASTEP"
    is_dft = True

    def __init__(self, settings=None, castepexe="castep.mpi"):
        super().__init__(settings)
        s = self.settings
        self.cascmd = s.exe or castep_command(nprocs=s.nprocs, ntasks=s.ntasks, nnodes=s.nnodes, mpirun=s.mpirun, castepexe=castepexe)

    def run(self, atoms, name="config", directory="."):
        s = self.settings
        E, atoms, P = call_castep(atoms, label=name, directory=os.path.join(directory, 'CASTEP-tmp'), cascmd=self.cascmd,
                                  inpseed=s.inpseed, base_dir=s.base_dir, pp_path=s.castep_pp_path)
        return float(E), atoms, P

    def describe(self):
        return "CASTEP (command: %s)" % self.cascmd


def get_engine(prog="ase", settings=None, **kwargs):
    """Create an engine: 'vasp', 'castep', 'mlip' (settings.mlip, e.g. 'mace') or 'ase' (settings.ase_pot: LJ, Morse, EMT)."""
    settings = settings or CalcSettings()
    key = prog.strip().lower()
    if key == "vasp":
        return VaspEngine(settings, **kwargs)
    if key == "castep":
        return CastepEngine(settings, **kwargs)
    if key == "mlip":
        if not settings.mlip:
            raise ValueError("No machine-learning potential was given (--mlip), e.g. --mlip mace")
        factory = lazy_calculator_factory(settings.mlip, model=settings.mlip_model, device=settings.device, **(settings.mlip_kwargs or {}))
        return AseEngine(factory, label=settings.mlip, settings=settings)
    if key in ("ase", "lj", "morse", "emt"):
        pot = settings.ase_pot if key == "ase" else key
        label = {"lj": "LJ", "morse": "Morse", "emt": "EMT"}.get(pot.lower(), pot)
        return AseEngine(lambda: classical_calculator(pot), label=label, settings=settings)
    raise ValueError("Unknown program '%s'; choose from vasp, castep, mlip, ase" % prog)
