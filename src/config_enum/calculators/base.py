"""Settings and the common interface of the calculation engines."""
from dataclasses import dataclass, field


@dataclass
class CalcSettings:
    """Calculation settings; the comments name the command-line option."""
    # VASP / CASTEP
    exe: str = None               # -vexe / -cexe (or $VASP_COMMAND / $CASTEP_COMMAND)
    idir: str = "INPUTS"          # -id: folder with INCAR, KPOINTS (and POTCAR) for VASP
    inpseed: str = None           # -iseed: seed of the .param/.cell files for CASTEP
    xc: str = "pbe"               # -xc
    magmoms: list = None          # -mgm, e.g. ['Fe', 5.0, 'Nb', 0.6]
    hubU: dict = field(default_factory=dict)  # -hubU, e.g. {'Fe': {'d': 2.7}}
    makepotcar: bool = False      # -mp
    potcarxc: str = "potpaw_PBE"  # --potcarxc
    nprocs: int = 32              # -np
    ntasks: int = 1               # -nt
    nnodes: int = 1               # -nn
    mpirun: bool = False          # -mpi
    castep_pp_path: str = None    # $CASTEP_PP_PATH (default: $HOME/APPS/pspots)
    base_dir: str = "."           # folder from which the run was started (where idir and the .usp/.param files are)
    # ASE calculators (classical and machine-learning potentials)
    ase_pot: str = "LJ"           # -pot
    optimise: bool = False        # -opt
    fmax: float = 0.05            # --fmax
    steps: int = 200              # --steps
    optimiser: str = "BFGS"       # --optimiser
    relax_cell: bool = False      # --relax_cell
    mlip: str = None              # --mlip
    mlip_model: str = None        # --mlip_model
    device: str = "auto"          # --device
    mlip_kwargs: dict = field(default_factory=dict)  # --mlip_kwargs


class Engine:
    """A calculation engine. run() returns (energy in eV, atoms, pressure in GPa)."""
    label = "E"
    is_dft = False

    def __init__(self, settings=None):
        self.settings = settings or CalcSettings()

    def run(self, atoms, name="config", directory="."):
        raise NotImplementedError

    def describe(self):
        return self.label
