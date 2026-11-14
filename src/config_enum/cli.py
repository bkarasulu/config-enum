"""Command line of config_enum (``config_enum`` command, or ``python config_enum.py``)."""
import argparse
import ast
import os
import sys

from . import __version__
from .calculators import CalcSettings, get_engine
from .calculators.classical import ASE_POTS
from .calculators.mlip import BACKENDS as MLIP_BACKENDS
from .calculators.optimise import OPTIMISERS
from .log import say, setup_logging
from .utils import parse_ranges
from .workflow import EnumWorkflow, WorkflowConfig

# Common variables.
inpdir = "INPUTS"
seed = 'config'
outdir = "OUT-FILES"


def _potcar_dirs():
    try:
        vasppp = os.environ['VASP_PP_PATH']
        ppdirs = os.listdir(vasppp)
        ppdirs.sort()
    except (KeyError, OSError):
        ppdirs = None
    return ppdirs or None


def parse_hubbard_u(hubU):
    """Parse the Hubbard_U info and make the dict: ['Fe', 'd', '2.7,', 'Sm', 'f', '6.1'] -> {'Fe': {'d': 2.7}, 'Sm': {'f': 6.1}}"""
    d = {}
    if hubU:
        U = " ".join(hubU) if not isinstance(hubU, str) else hubU
        for xx in U.split(','):
            x = xx.split()
            if len(x) == 0:
                continue
            if len(x) != 3:
                raise ValueError("Cannot parse -hubU '%s'; expected e.g. -hubU Fe d 2.7, Sm f 6.1" % U)
            d.setdefault(x[0], {})[x[1]] = float(x[2])
    return d


def parse_no_dopants(noDopants):
    """Get the dopant/vacancy number list: ['0:4', '24'] -> [0, 1, 2, 3, 4, 24]"""
    try:
        nDops = sorted(set(parse_ranges([str(y) for y in noDopants])))
    except ValueError:
        raise ValueError("Cannot parse -nd '%s'; expected e.g. -nd 0:4 24" % " ".join(str(y) for y in noDopants))
    if len(nDops) == 0 or min(nDops) < 0:
        raise ValueError("-nd: the number of dopants/vacancies must be given (0 or larger)")
    return nDops


def build_parser():
    parser = argparse.ArgumentParser(
        prog="config_enum",
        description='This script generates the input structures required for a configuration enumeration (including vacancies and dopants) '
                    'for any software supported by ASE, and computes their energies with VASP, CASTEP, machine-learning potentials (MLIPs) or '
                    'simple classical potentials. \nRequires ASE and SPGlib packages.\n Sample usage: config_enum.py -i LLZO-Ia-3dCollCode261302.cif '
                    '-prim -at Li -wl 96h  -dt Va -nd 0:2  -sf -ot vasp -o POSCAR')

    parser.add_argument('-i', '--inpf', type=str, required=True, nargs='*', help='Input file(s), any file type compatible with ASE.')
    parser.add_argument('-o', '--outf', type=str, required=False, help='Output file name. Def: Automatic naming: e.g. %s_ID.' % seed)
    parser.add_argument('-it', '--itype', type=str, required=False, help='Input file type. Def: determined automatically from the extension.')
    parser.add_argument('-ot', '--otype', default='res', type=str, required=False, help='Output file type, default: .res (SHELX format)')
    parser.add_argument('-od', '--odir', default=outdir, type=str, required=False, help='Output directory (kept for compatibility; not used: the output is written to StepX/input_N folders).')

    parser.add_argument('-ow', '--overwrite', default=False, action='store_true', help='overwrite if output folder exists. Def: No')
    parser.add_argument('-rs', '--restart', default=False, action='store_true', help='Do a continuation run (all converged VASP calculations will be read in and will not be redone). Def: No')
    parser.add_argument("-iseed", "--inpseed", type=str, default=None, help="Seed name for reading the CASTEP .param and .cell files in order to use the settings within.")
    parser.add_argument('-sf', '--sepfol', default=False, action='store_true', help='to save geometries in separate folders (e.g. for VASP runs). Default: all output of an input structure is written in the same directory.')

    parser.add_argument('-prim', '--ifPrim', default=False, action='store_true', help='use the primitive cell instead of supercell (requires SPGlib installed). Highly recommended for computational efficiency. Def: false')
    parser.add_argument('-ref', '--initCrysFile', type=str, help='File for obtaining the initial crystal site (symmetry) information. Can be any type that ASE supports; it must have the same cell and atom order as the structure to be doped. Def: same as the input file (-i)')
    parser.add_argument('-tol', '--tol', type=float, default=1e-3, help="The symmetry tolerance. Def: 1e-3")
    parser.add_argument('-nosymm', '--nosymm', default=False, action='store_true', help='switch off the use of symmetry info, e.g. equivalent sites (not recommended).')
    parser.add_argument('-r', '--rep', nargs=3, type=int, help="Repeat the input structure in a, b, c directions to get a supercell, e.g. -r 2 2 1")

    parser.add_argument('-nl', '--list', help='atom list(count starts at 1); e.g. -nl 1,2,5,8:12')
    parser.add_argument('-at', '--aType', nargs='*', type=str, help='atom type to be included in enumeration; e.g. -at Na Li')
    parser.add_argument('-wl', '--wycklist', nargs='*', type=str, help='List of Wyckoff positions of crys. sites to include in enumeration. A list of available sites can be obtained via dry run (-dry); e.g. -wl 24d 12a')
    parser.add_argument('-excl', '--excludelist', nargs='*', type=str, help='List of crystal sites (atom indices) to exclude from the configuration enumeration (these atoms will not be considered even though selected by -nl, -at or -wl earlier).  e.g. 1:4 24  for [1,2,3,4,24], atom numbering starts at 1.')
    parser.add_argument('-dexcl', '--dopeexcluded', default=False, action='store_true', help='Add dopant/vacancy to the atoms/crys-sites given in the excluded list (in all configurations). This option is good for doing iterative vacancy/dopant addition with the previous knowledge of lower dopant/vacancy compositions. Def: False.')
    parser.add_argument('-dt', '--dType', required=True, type=str, help='Dopant type to consider in the enumeration; This can be an atom symbol (e.g. Na) or vacancy (e.g. Va)')
    parser.add_argument('-dexcl_type', '--dexcl_type', required=False, type=str, help='Dopant type to consider for the sites given in the excluded list (see -dexcl keyword for details). Def: equal to dType', default="")
    parser.add_argument('-nd', '--noDopants', required=True, nargs='*', type=str, help='(List of) number of dopant/vacancy to consider in enumeration; e.g. 0:4 24  for 0,1,2,3,4 and 24 vacancies)')

    parser.add_argument('-mod', '--modify', type=str, help="To make some final touches on the configurations selected for enumeration. You can combine multiple modifications using multiple statements separated by ';', the order of the statements will be strictly followed. This is a useful option to do necessary modifications to compensate for the dopant effect  e.g. balancing the total system charge  by deleting additional (nearest) cations upon introducing higher-valence cation(s). The reference can be an atom type or 'dopant' (the doped/vacant sites). Def: No additional modifications.  Example usage: -mod 'replace all within 2 of Al with Va' or -mod 'remove Li nearest 1 of dopant'.")
    parser.add_argument('-ai', '--addint', default=False, action='store_true', help='To add the interstices as X atoms (select them with -at X) based on the Wyckoff positions data for the reference symmetry group. Default: False.')
    parser.add_argument('-ase', '--ASE_symm', default=False, action='store_true', help='Do a final check of the symmetrically-unique configurations by comparing the whole structures with ASE (SymmetryEquivalenceCheck). Default: False.')

    parser.add_argument('-castep', '--castep', default=False, action='store_true', help='Run quick, low-accuracy CASTEP calculations to estimate the configuration energies. Default: a simple ASE potential (see -pot) is used.')
    parser.add_argument('-vasp', '--vasp', default=False, action='store_true', help='Use VASP to compute configuration energies. Inputs (INCAR, KPOINTS, POTCAR) are read from the input directory (-id)')
    parser.add_argument('-vexe', '--vaspexe', type=str, required=False, default='vasp_std', help='Vasp executable. Def: vasp_std')
    parser.add_argument('-cexe', '--castepexe', type=str, required=False, default='castep.mpi', help='CASTEP executable. Def: castep.mpi')

    parser.add_argument('-copy', '--copyinp', default=False, action='store_true', help='to copy the files contained in "%s" folder (with -sf). Default: False.' % inpdir)
    parser.add_argument('-id', '--idir', default=inpdir, type=str, required=False, help='Input directory to copy the files from, default: %s' % inpdir)

    parser.add_argument('-mp', '-makepotcar', '--makepotcar', default=False, action='store_true', help='to compile POTCAR (for VASP) using actual atomic content. The environment variable ($VASP_PP_PATH) should be defined beforehand. Default: False.')
    parser.add_argument('--potcarxc', default="potpaw_PBE", type=str, choices=_potcar_dirs(), help='XC functional to use in making of POTCAR. Def: potpaw_PBE')
    parser.add_argument('-xc', '--xc', type=str, required=False, default='pbe', help='Exchange correlation functional to use (pseudopots will be selected accordingly, e.g. LDA, PBE,PBEsol, etc. Def: PBE')
    parser.add_argument("-hubU", "--hubU", type=str, nargs="*", help="For defining the Hubbard U parameters (in eV) for specific atom and orbital types, e.g. -hubU Fe d 2.7, Sm f 6.1")
    parser.add_argument("-mgm", "--magmoms", default=None, nargs="*", required=False, help="Magnetic moments for a collinear calculation. Eg, 'Fe 5.0 Nb 0.6 O 0.6' Def.: None. If a defined element is not present in the POSCAR, no MAGMOM will be set for it.")

    parser.add_argument('-np', '--nprocs', type=int, default=32, help="No of processes to start for each CASTEP/VASP calculation through srun/mpirun. Def:32")
    parser.add_argument('-nt', '--ntasks', type=int, default=1, help="No of calculations (tasks) to run simultaneously. Def:1")
    parser.add_argument('-nn', '--nnodes', type=int, default=1, help="No of nodes to run CASTEP/VASP runs through srun. Def:1")
    parser.add_argument('-mpi', '--mpirun', default=False, action='store_true', help='Use mpirun for parallel runs. Default: srun is used')

    parser.add_argument("-t", "--temp", type=float, default=298.15, help="Temperature at which the configurational thermodynamics is computed. Def: 298.15K")

    parser.add_argument("-seq", "--seq", type=int, default=0, help="If to run an automated sequential doping to reach the target number of dopants (-nd): the number of dopants/vacancies added at each step. The configurations are selected at each step based on energy/symmetry (--seqtype). Def: 0 (no sequential doping; all dopants are added at once).")
    parser.add_argument("-spick", "--seqpick", type=int, default=1, help="How many configs to pick at each step in a sequential doping. Def: 1")
    parser.add_argument("-stype", "--seqtype", type=str, choices=['E', 'S', 'D', 'R'], default='E', help="How to sort the configurations at each step of  a sequential doping. Options: 'E': Based on Energy (low-energy configs are picked); 'S': based on symmetry (highest-symmetry configs picked); 'D': config(s) with highest degeneracy are picked; 'R': configs are picked randomly. Def: 'E' ")

    parser.add_argument('-v', '--verb', default=False, action='store_true', help='Verbose output. Default: False.')
    parser.add_argument('-dry', '--dry', default=False, action='store_true', help='Dry run to get the equivalent (Wyckoff) positions. Default: False.')
    parser.add_argument('-noskip', '--noskip', default=False, action='store_true', help='Do not skip problematic configs (and assigning them with 0 energy and P1 symmetry). Default: Bad configs are skipped.')

    # Options added in the refactored version
    parser.add_argument("-ssym", "--seqsymm", type=str, choices=['current', 'ref'], default='current', help="Symmetry operations used after the first step of a sequential doping: 'current': those of the current (doped) structure; 'ref': those of the Step 1 (reference) structure, as in the original script (fewer configurations, but inequivalent ones can be merged). Def: current")
    parser.add_argument("--wyckoff_file", type=str, default=None, help="File with Wyckoff positions (in the format of the wyckoff.dat of the original script) to use for -ai instead of the table of the package. Def: $CONFIG_ENUM_WYCKOFF, or the table of the package")
    parser.add_argument('-q', '--quiet', default=False, action='store_true', help='Only show warnings and errors. Default: False.')
    parser.add_argument("--version", action="version", version="config_enum %s" % __version__)

    # ASE calculators (classical and machine-learning potentials)
    parser.add_argument("-pot", "--ase_pot", default='LJ', type=str, choices=list(ASE_POTS), help="Simple potential of ASE used if none of -vasp, -castep and --mlip is given (for quick tests only; LJ and Morse do not distinguish the atom types). Def: LJ")
    parser.add_argument("-opt", "--optimise", default=False, action='store_true', help="ASE/MLIP: optimise the geometry of each configuration (Def: single-point energies). The structures passed on to the next step of a sequential doping are always the unoptimised ones.")
    parser.add_argument("--fmax", default=0.05, type=float, help="ASE/MLIP: force convergence criterion of the geometry optimisation. Def: 0.05 eV/A")
    parser.add_argument("--steps", default=200, type=int, help="ASE/MLIP: max number of geometry optimisation steps. Def: 200")
    parser.add_argument("--optimiser", default="BFGS", type=str.upper, choices=list(OPTIMISERS), help="ASE/MLIP: ASE optimiser. Def: BFGS")
    parser.add_argument("--relax_cell", default=False, action='store_true', help="ASE/MLIP: also relax the cell in the geometry optimisation (implies -opt). Def: False")

    # Machine-learning potentials
    parser.add_argument("-mlip", "--mlip", default=None, type=str, help="Use a machine-learning potential to compute the configuration energies: %s, or custom:package.module:function (any function returning an ASE calculator). Def: none" % ", ".join(MLIP_BACKENDS))
    parser.add_argument("--mlip_model", default=None, type=str, help="Model name or path for --mlip, e.g. 'medium-mpa-0' or 'my.model' (MACE), '7net-0' (SevenNet), 'uma-s-1' (fairchem). Def: the backend default")
    parser.add_argument("--device", default="auto", type=str, help="Device for the machine-learning potential: cpu, cuda, mps or auto. Def: auto (GPU if available)")
    parser.add_argument("--mlip_kwargs", default=None, type=str, help="Extra keyword arguments for the ML potential calculator as a Python dict, e.g. \"{'default_dtype': 'float32'}\" or \"{'task_name': 'omat'}\". Def: none")

    # Input files
    parser.add_argument("--input", default=None, type=str, help="Input file (TOML, YAML or JSON) with options; keys are the long option names (e.g. inpf, aType, dType, noDopants). Options given on the command line override the file.")
    parser.add_argument("--save_input", default=None, type=str, help="Write the effective options of this run to a TOML/YAML/JSON file (by extension) that can be used with --input.")
    return parser


INPUT_EXCLUDE = {"help", "version", "input", "save_input"}


def load_input_file(path):
    """Read an input file (TOML, YAML or JSON, by extension) into a dict."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError:
            raise ValueError("YAML input files need PyYAML (pip install pyyaml); TOML and JSON work without it.")
        with open(path) as f:
            data = yaml.safe_load(f) or {}
    elif ext == ".json":
        import json
        with open(path) as f:
            data = json.load(f)
    else:
        try:
            import tomllib
        except ImportError:
            import tomli as tomllib  # Python < 3.11
        with open(path, "rb") as f:
            data = tomllib.load(f)
    if not isinstance(data, dict):
        raise ValueError("Input file %s must contain a table/mapping of options." % path)
    # allow optional sections, e.g. [structures], [calculation] -> flatten
    flat = {}
    for k, v in data.items():
        if isinstance(v, dict) and k not in ("mlip_kwargs",):
            flat.update(v)
        else:
            flat[k] = v
    return flat


def _option_index(parser):
    # name (dest, long or short option without dashes) -> argparse action
    idx = {}
    for a in parser._actions:
        if a.dest in INPUT_EXCLUDE or a.dest == argparse.SUPPRESS:
            continue
        idx.setdefault(a.dest, a)
        for opt in a.option_strings:
            idx.setdefault(opt.lstrip("-"), a)
    return idx


def _convert(action, value):
    """Convert a value from an input file to what argparse would produce for this option."""
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        if not isinstance(value, bool):
            raise ValueError("option '%s' expects true or false, got %r" % (action.option_strings[-1].lstrip("-"), value))
        return value if isinstance(action, argparse._StoreTrueAction) else (not value)

    def one(v):
        if isinstance(v, dict):
            v = repr(v)  # e.g. mlip_kwargs (parsed later with ast.literal_eval)
        if isinstance(v, bool):
            raise ValueError("option '%s' does not take true/false" % action.dest)
        if action.type is None:
            return v
        return action.type(v)
    if action.nargs in ("+", "*") or isinstance(action.nargs, int):
        if isinstance(value, str):
            value = value.split()
        elif not isinstance(value, (list, tuple)):
            value = [value]
        return [one(v) for v in value]
    if isinstance(value, (list, tuple)):
        value = ",".join(str(v) for v in value)  # e.g. list = [1, 2, 5]
    return one(value)


def apply_input_file(parser, path):
    """Use the options of an input file as defaults of the parser (command-line options still override them)."""
    values = load_input_file(path)
    idx = _option_index(parser)
    unknown = sorted(k for k in values if k not in idx)
    if unknown:
        raise ValueError("Unknown option(s) in %s: %s" % (path, ", ".join(unknown)))
    defaults = {}
    for k, v in values.items():
        a = idx[k]
        defaults[a.dest] = _convert(a, v)
        a.required = False  # given in the file
    parser.set_defaults(**defaults)
    return defaults


def _toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        return '"%s"' % v.replace("\\", "\\\\").replace('"', '\\"')
    if isinstance(v, (list, tuple)):
        return "[%s]" % ", ".join(_toml_value(x) for x in v)
    raise ValueError("Cannot write %r to TOML" % (v,))


def save_input_file(args, path, parser=None):
    """Write the parsed options (those differing from the defaults, plus the required ones) to a TOML/YAML/JSON file."""
    parser = parser or build_parser()
    defaults = {a.dest: parser.get_default(a.dest) for a in parser._actions}
    data = {}
    for k, v in sorted(vars(args).items()):
        if k in INPUT_EXCLUDE or v is None:
            continue
        if k in ("inpf", "dType", "noDopants") or v != defaults.get(k):
            data[k] = list(v) if isinstance(v, tuple) else v
    ext = os.path.splitext(path)[1].lower()
    with open(path, "w") as f:
        if ext == ".json":
            import json
            json.dump(data, f, indent=1)
        elif ext in (".yaml", ".yml"):
            import yaml
            yaml.safe_dump(data, f, sort_keys=True)
        else:
            f.write("# config_enum input file (use with: config_enum --input %s)\n" % os.path.basename(path))
            for k, v in data.items():
                f.write("%s = %s\n" % (k, _toml_value(v)))
    return data


def get_prog(args):
    """Which engine is used: vasp, castep, mlip or ase."""
    if args.vasp:
        return 'vasp'
    if args.castep:
        return 'castep'
    if args.mlip:
        return 'mlip'
    return 'ase'


def make_settings(args, base_dir='.'):
    """CalcSettings from the parsed command-line arguments."""
    mlip_kwargs = {}
    if args.mlip_kwargs:
        try:
            mlip_kwargs = ast.literal_eval(args.mlip_kwargs)
        except (ValueError, SyntaxError):
            raise ValueError("Cannot parse --mlip_kwargs '%s' (expected a Python dict)" % args.mlip_kwargs)
        if not isinstance(mlip_kwargs, dict):
            raise ValueError("--mlip_kwargs must be a Python dict, got '%s'" % args.mlip_kwargs)
    if args.magmoms and len(args.magmoms) % 2:
        raise ValueError("-mgm needs pairs of atom type and magnetic moment, e.g. -mgm Fe 5.0 O 0.6")
    return CalcSettings(idir=args.idir, inpseed=args.inpseed, xc=args.xc, magmoms=args.magmoms, hubU=parse_hubbard_u(args.hubU),
                        makepotcar=args.makepotcar, potcarxc=args.potcarxc, nprocs=args.nprocs, ntasks=args.ntasks, nnodes=args.nnodes,
                        mpirun=args.mpirun, base_dir=os.path.abspath(base_dir), ase_pot=args.ase_pot,
                        optimise=args.optimise or args.relax_cell, fmax=args.fmax, steps=args.steps, optimiser=args.optimiser,
                        relax_cell=args.relax_cell, mlip=args.mlip, mlip_model=args.mlip_model, device=args.device, mlip_kwargs=mlip_kwargs)


def config_from_args(args, base_dir='.', command=''):
    """WorkflowConfig from the parsed command-line arguments."""
    return WorkflowConfig(
        dType=args.dType, noDopants=parse_no_dopants(args.noDopants), outf=args.outf, itype=args.itype, otype=args.otype,
        overwrite=args.overwrite, restart=args.restart, sepfol=args.sepfol, ifPrim=args.ifPrim, initCrysFile=args.initCrysFile, tol=args.tol,
        nosymm=args.nosymm, rep=tuple(args.rep) if args.rep else None, list=args.list, aType=args.aType, wycklist=args.wycklist,
        excludelist=args.excludelist, dopeexcluded=args.dopeexcluded, dexcl_type=args.dexcl_type, modify=args.modify, addint=args.addint,
        wyckoff_file=args.wyckoff_file, ASE_symm=args.ASE_symm, copyinp=args.copyinp, idir=args.idir,
        makepotcar=args.makepotcar, potcarxc=args.potcarxc, ntasks=args.ntasks, temp=args.temp, seq=args.seq, seqpick=args.seqpick,
        seqtype=args.seqtype, seqsymm=args.seqsymm, verb=args.verb, dry=args.dry, noskip=args.noskip, base_dir=os.path.abspath(base_dir),
        command=command)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()

    # first pass: only look for --input, to use the file as defaults
    pre = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    pre.add_argument("--input", default=None)
    known, _ = pre.parse_known_args(argv)
    if known.input:
        try:
            apply_input_file(parser, known.input)
        except (OSError, ValueError) as err:
            parser.error(str(err))
    args = parser.parse_args(argv)
    setup_logging(-1 if args.quiet else (1 if args.verb else 0))
    say()

    if args.save_input:
        save_input_file(args, args.save_input, parser)
        say("Options written to %s" % args.save_input)
    if not args.inpf:
        parser.error("no input file was given (-i)")
    if not (args.aType or args.wycklist or args.list) and not args.dry:
        parser.error("The list of atoms (and their corresponding crystalline sites) which will be replaced by vacancy or dopants must be "
                     "given. Use -nl, -at or -wl for this purpose.")

    try:
        settings = make_settings(args)
        config = config_from_args(args, command="config_enum " + " ".join(argv))
        if config.seq > 0 and len(config.noDopants) > 1:
            raise ValueError('Multiple no of dopant entries (-nd) is not supported with -seq option...')
    except (ValueError, SyntaxError) as err:
        parser.error(str(err))
    if settings.hubU:
        say("Hubbard U values to be used: ", settings.hubU)

    prog = get_prog(args)
    try:
        if prog == 'vasp':
            engine = get_engine(prog, settings, vaspexe=args.vaspexe)
        elif prog == 'castep':
            engine = get_engine(prog, settings, castepexe=args.castepexe)
        else:
            engine = get_engine(prog, settings)
        if prog == 'mlip' and not args.dry:  # load the model now, so that a missing package or model is reported before any work is done
            engine.get_calculator()
    except (ImportError, ValueError) as err:
        parser.error(str(err))
    if prog in ('vasp', 'castep'):
        say("Total of %d cores will be used in %d parallel tasks." % (args.nprocs, args.ntasks))
        say(engine.describe())

    try:
        EnumWorkflow(engine, config).run(args.inpf)
    except (RuntimeError, FileNotFoundError, ValueError) as err:
        say("%s\nTerminating..." % err)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
