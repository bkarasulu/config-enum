"""The configuration-enumeration workflow: (sequential) doping, energies and the selection of configurations.

For every step (one step unless sequential doping, -seq, is used) and every input structure:
    1. prepare the structure (interstices, primitive cell, supercell) and its symmetry operations,
    2. select the atoms/crystal sites and enumerate the symmetrically-unique configurations,
    3. write the structures,
    4. compute the energies (VASP, CASTEP, a machine-learning potential or a classical potential),
    5. configurational thermodynamics (degeneracies, reduced energies and probabilities),
and finally the configurations are sorted and picked for the next step.
"""
import glob
import json
import multiprocessing
import os
import re
import shutil
import time
from dataclasses import dataclass, field

import numpy as np
from ase.build import sort

from .calculators.vasp import make_potcar
from .enumeration import (PARENT_ARRAY, apply_dopants, canonical_keys, changed_sites, enumerate_configs, parse_exclude_list,
                          parse_modify, select_atoms, set_parent_sites)
from .io import partial_occupancies, read_struct, write_res_result, write_struct
from .log import debug, say, warn
from .symmetry import (find_interstices, find_primitive_cell, get_equivalent_sites, get_symmetry, identity_symmetry,
                       print_equivalent_sites, spacegroup_number, spacegroup_string)
from .thermo import boltz_dist
from .utils import hashid

seed = 'config'  # common name of the configurations


@dataclass
class WorkflowConfig:
    """Options of the workflow; the comments name the command-line option."""
    dType: str = 'Va'              # -dt
    noDopants: list = field(default_factory=lambda: [1])  # -nd (list of numbers)
    outf: str = None               # -o
    itype: str = None              # -it
    otype: str = 'res'             # -ot
    overwrite: bool = False        # -ow
    restart: bool = False          # -rs
    sepfol: bool = False           # -sf
    ifPrim: bool = False           # -prim
    initCrysFile: str = None       # -ref
    tol: float = 1e-3              # -tol
    nosymm: bool = False           # -nosymm
    rep: tuple = None              # -r
    list: str = None               # -nl
    aType: list = None             # -at
    wycklist: list = None          # -wl
    excludelist: list = None       # -excl
    dopeexcluded: bool = False     # -dexcl
    dexcl_type: str = ''           # -dexcl_type
    modify: str = None             # -mod
    addint: bool = False           # -ai
    wyckoff_file: str = None       # --wyckoff_file
    ASE_symm: bool = False         # -ase
    copyinp: bool = False          # -copy
    idir: str = 'INPUTS'           # -id
    makepotcar: bool = False       # -mp
    potcarxc: str = 'potpaw_PBE'   # --potcarxc
    ntasks: int = 1                # -nt
    temp: float = 298.15           # -t
    seq: int = 0                   # -seq
    seqpick: int = 1               # -spick
    seqtype: str = 'E'             # -stype
    seqsymm: str = 'current'       # -ssym
    verb: bool = False             # -v
    dry: bool = False              # -dry
    noskip: bool = False           # -noskip
    base_dir: str = '.'            # folder of the run (the StepX folders are created here)
    command: str = ''              # the command line (written to the 'degeneracies' file)


@dataclass
class InputStructure:
    """An input structure of a step: a file, or a configuration picked in the previous step."""
    path: str                      # file name (as shown in the output)
    atoms: object = None           # the structure (read from path if None)
    atoms_list: list = None        # selected sites still available (None: determined from the options)
    root: int = None               # index of the Step 1 input structure it derives from


@dataclass
class ConfigResult:
    """A symmetrically-unique configuration."""
    name: str = ''                 # config_<ID>
    step: int = 1
    input_index: int = 0           # index of the input structure of the step
    nD: int = 0                    # number of dopants/vacancies added in this step
    sites: list = field(default_factory=list)  # indices of the doped atoms in the input structure
    degeneracy: int = 1
    atoms: object = None           # initial structure (dopants on the lattice sites, not optimised)
    atoms_list: list = field(default_factory=list)  # selected sites still available (indices in next_atoms)
    next_atoms: object = None      # the structure passed on to the next step (as atoms, plus the unused interstices)
    file: str = ''                 # structure file, relative to the run folder
    energy: float = None
    E_reduced: float = None        # energy - kT ln(degeneracy)
    prob: float = None
    pressure: float = 0.
    spacegroup: str = 'P1 (1)'
    final_atoms: object = None     # structure after the calculation
    result_file: str = ''          # .res file with the final structure and energy
    failed: bool = False
    picked: bool = False
    root: int = 0                  # index of the Step 1 input structure it derives from
    key: tuple = None              # the same for symmetrically-equivalent structures (see enumeration.canonical_keys)

    def as_dict(self):
        d = {k: getattr(self, k) for k in ('name', 'step', 'input_index', 'nD', 'degeneracy', 'file', 'energy', 'E_reduced',
                                           'prob', 'pressure', 'spacegroup', 'result_file', 'failed', 'picked')}
        d['sites'] = [int(i) for i in self.sites]
        d['formula'] = self.atoms.get_chemical_formula('hill') if self.atoms is not None else ''
        return d


@dataclass
class StepResult:
    step: int = 1
    nDops: list = field(default_factory=list)
    inputs: list = field(default_factory=list)     # file names
    configs: list = field(default_factory=list)    # all ConfigResult of the step
    picked: list = field(default_factory=list)     # the ConfigResult picked for the next step
    summary: str = ''


@dataclass
class WorkflowResults:
    steps: list = field(default_factory=list)

    @property
    def configs(self):
        return [c for s in self.steps for c in s.configs]

    @property
    def picked(self):
        return self.steps[-1].picked if self.steps else []

    def as_dict(self):
        return {'steps': [{'step': s.step, 'nDops': s.nDops, 'inputs': s.inputs, 'configs': [c.as_dict() for c in s.configs],
                           'picked': [c.file for c in s.picked]} for s in self.steps]}


# ----------------------------------------------------------------------------------------- parallel calculations
_WORK = {}


def _run_one(i):
    """Run the calculation of configuration i (of the list in _WORK); returns (i, energy, atoms, pressure, error)."""
    engine = _WORK['engine']
    conf = _WORK['configs'][i]
    atoms = conf.atoms.copy()
    try:
        E, atoms, P = engine.run(atoms, name=conf.name, directory=_WORK['directory'])
        atoms.calc = None
        return i, E, atoms, P, None
    except (KeyboardInterrupt, ImportError):
        raise
    except Exception as err:
        atoms.calc = None
        return i, 0.0, atoms, 0., "%s: %s" % (type(err).__name__, err)


def run_calculations(engine, configs, directory, ntasks=1):
    """Run the engine on all configurations (ntasks at a time); yields (i, energy, atoms, pressure, error) in order."""
    _WORK.update(engine=engine, configs=configs, directory=directory)
    try:
        if ntasks > 1 and len(configs) > 1:
            with multiprocessing.get_context("fork").Pool(ntasks) as pool:
                yield from pool.imap(_run_one, range(len(configs)))
        else:
            for i in range(len(configs)):
                yield _run_one(i)
    finally:
        _WORK.clear()


def ase_symmetry_filter(configs):
    """A last check of the pre-filtered set of configs, comparing the whole structures with ASE (this also considers
    different atomic orders and all atoms in the comparison). The degeneracies of identical configs are added."""
    from ase.utils.structure_comparator import SymmetryEquivalenceCheck as ASE_SEC
    say('compare structures using ASE')
    comp = ASE_SEC()
    uniq_configs = []
    found = []
    for i in range(len(configs)):
        if i in found:
            continue
        found.append(i)
        for j in range(i + 1, len(configs)):
            if j in found or configs[i].nD != configs[j].nD:
                continue
            if comp.compare(configs[i].atoms, configs[j].atoms):  # 1st ref, 2nd compared.
                say('%d matches to %d ' % (i, j))
                found.append(j)
                configs[i].degeneracy += configs[j].degeneracy
        uniq_configs.append(i)
    say('Unique list from ASE: ', uniq_configs)
    return [configs[i] for i in uniq_configs]


# ------------------------------------------------------------------------------------------------------ workflow
class EnumWorkflow:
    """The enumeration workflow.

    Example::

        from config_enum.calculators import CalcSettings, get_engine
        from config_enum.workflow import EnumWorkflow, WorkflowConfig

        engine = get_engine("mlip", CalcSettings(mlip="mace", mlip_model="small", optimise=True))
        cfg = WorkflowConfig(aType=["Li"], dType="Va", noDopants=[2])
        results = EnumWorkflow(engine, cfg).run(["LiCoO2.cif"])
    """

    def __init__(self, engine, config):
        self.engine = engine
        self.cfg = config
        self.base = os.path.abspath(config.base_dir)
        self.results = WorkflowResults()
        self._ref_symmetry = {}  # (all_W, all_w) of step 1 for each input file
        self._parents = {}       # the Step 1 structures (before doping), by input file

    # ------------------------------------------------------------------------------------------------- helpers
    def step_dir(self, step):
        return os.path.join(self.base, 'Step%d' % step)

    def _check_existing(self):
        cfg = self.cfg
        if os.path.exists(self.step_dir(1)):
            if cfg.overwrite:
                say("StepX folders exist. Overwriting as requested.")
                for fn in os.listdir(self.base):
                    if re.fullmatch(r'Step\d+', fn) and os.path.isdir(os.path.join(self.base, fn)):
                        shutil.rmtree(os.path.join(self.base, fn))
            elif cfg.restart:
                say("StepX folders exist. Will try to restart the calculations as requested.")
            else:
                raise RuntimeError("StepX folders exists. Not overwriting, use -ow to overwrite (or -rs to restart).")

    def get_reference(self, atoms):
        """The structure used for the symmetry information (-ref), or the structure itself."""
        cfg = self.cfg
        if cfg.initCrysFile:  # if an input file for crys. site info is provided.
            # Use an external ref. for symm ops. to determine the Wyckoff pos.
            ref = cfg.initCrysFile if os.path.isabs(cfg.initCrysFile) else os.path.join(self.base, cfg.initCrysFile)
            say('%s is used as reference for the symmetry information.' % cfg.initCrysFile)
            atoms_crys = read_struct(ref)
            if len(atoms_crys) != len(atoms) or not np.allclose(atoms_crys.cell.cellpar(), atoms.cell.cellpar(), atol=1e-2):
                warn("WARNING: the reference structure (%d atoms) and the structure to be doped (%d atoms) have different "
                     "cells; the symmetry operations and Wyckoff sites may not apply." % (len(atoms_crys), len(atoms)))
            return atoms_crys
        return atoms.copy()  # use input file as the template for symmetry.

    # ------------------------------------------------------------------------------------------------- the steps
    def prepare_structure(self, inp, step, find):
        """Read an input structure and apply the requested changes (interstices, primitive cell, supercell).

        Returns (atoms, rpt_no): rpt_no is the number of primitive cells in the original cell (with -prim).
        """
        cfg = self.cfg
        if inp.atoms is not None:
            atoms = inp.atoms.copy()
        else:
            atoms = read_struct(inp.path if os.path.isabs(inp.path) else os.path.join(self.base, inp.path), cfg.itype)
            partial = partial_occupancies(atoms)
            if partial:
                say("Partially occupied sites in %s (all sites are treated as occupied; use vacancies, -dt Va, to remove atoms):" % inp.path)
                for sym, o, n in partial:
                    say("  %s: %d sites with occupancy %.3g (%.1f atoms; %.1f to remove)" % (sym, n, o, n * o, n * (1 - o)))

        # Find the potential interstices/vacancies, checking the Wyckoff site list for the space group. This must be
        # done before finding the primitive and repetitions
        if step == 1 and cfg.aType and ('X' in cfg.aType and cfg.addint):
            atoms_crys = self.get_reference(atoms)
            sg = spacegroup_string(atoms_crys, cfg.tol)
            say(sg)
            atoms = find_interstices(atoms, tol=cfg.tol, verb=cfg.verb, sg=spacegroup_number(sg), wyckoff_file=cfg.wyckoff_file)

        rpt_no = 1
        if step == 1 and cfg.ifPrim:
            say("Finding the primitive cell using SPGlib as requested.")
            say("Space group= %s with tol=%.2e" % (spacegroup_string(atoms, cfg.tol), cfg.tol))
            atoms, rpt_no = find_primitive_cell(atoms, cfg.tol)

        if step == 1 and cfg.rep:
            atoms = atoms.repeat(cfg.rep)
            atoms = sort(atoms)
        return atoms, rpt_no

    def get_symmetry(self, atoms, step, root, inpf='', rpt_no=1):
        """Symmetry operations and the Wyckoff labels of the atoms: (all_W, all_w, wyckoffs)."""
        cfg = self.cfg
        if cfg.nosymm:  # The nosymm case
            all_W, all_w = identity_symmetry()
            return all_W, all_w, {}

        if step > 1 and cfg.seqsymm == 'ref' and root in self._ref_symmetry:
            all_W, all_w = self._ref_symmetry[root]
            say("%d symmetry operations of the reference (Step 1) structure are used." % len(all_W))
            return all_W, all_w, {}

        # In the later steps of a sequential doping, the symmetry of the current (doped) structure is used: only the
        # operations that leave the dopants/vacancies already present unchanged relate equivalent configurations.
        atoms_crys = self.get_reference(atoms) if step == 1 else atoms.copy()
        data = get_symmetry(atoms_crys, cfg.tol)
        if data is not None:
            all_W = data['rotations']
            all_w = data['translations']
            say("%d symmetry operations were detected. " % (len(all_w)))
            say("Space group= %s with tol=%.2e" % (data['spacegroup'], cfg.tol))
            debug(data['equivalent_atoms'])
        else:
            all_W, all_w = identity_symmetry()

        eSites, wyckoffs = get_equivalent_sites(atoms, data, rpt_no=rpt_no)
        print_equivalent_sites(eSites, inpf)
        return all_W, all_w, wyckoffs

    def select(self, atoms, wyckoffs, inp, step):
        """The atoms selected for the enumeration and the excluded list: (atoms_list, eList)."""
        cfg = self.cfg
        if inp.atoms_list is not None:  # sequential doping: the sites left from the previous step
            atoms_list = list(inp.atoms_list)
        else:
            # Parse the user input as to atom selection.
            atoms_list = select_atoms(atoms, aType=cfg.aType, wycklist=cfg.wycklist, wyckoffs=wyckoffs, atom_list=cfg.list)
        if len(atoms_list) == 0:
            raise RuntimeError("No atoms were selected based on your input. Check the atom selection options (i.e. -nl, -at and -wl). "
                               "See help (-h) for further information.")

        say("\n%d atoms selected." % len(atoms_list))
        say("List of atoms (index starts at 0) selected for the enumeration: ", atoms_list)

        eList = []
        if cfg.excludelist and step == 1:
            eList = [x for x in parse_exclude_list(cfg.excludelist) if x < len(atoms)]
            say("Excluded list: ", eList)
            atoms_list = [x for x in atoms_list if x not in eList]
            say("New list of selected %d atoms: " % len(atoms_list), atoms_list)
            if cfg.dopeexcluded:
                if len(eList) == 0:
                    say("Exclude list is empty, switching off -dexcl keyword.")
                else:
                    say("Atoms in the excluded list will be doped/left vacant (%s) in all configurations as requested."
                        % (cfg.dexcl_type or cfg.dType))
        return atoms_list, eList

    def generate(self, atoms, atoms_list, eList, nDops, all_W, all_w, step, find):
        """Create the symmetrically-unique configurations (ConfigResult objects, not yet written)."""
        cfg = self.cfg
        mods = parse_modify(cfg.modify, atomSyms=list(np.unique(atoms.get_chemical_symbols())))
        dexcl = eList if cfg.dopeexcluded else []
        configs = []
        for nD in nDops:
            if nD > len(atoms_list):
                say("%d vacancies/dopants requested, but only %d atoms are selected, skipping..." % (nD, len(atoms_list)))
                continue
            initT1 = time.time()
            confs, degens = enumerate_configs(atoms, atoms_list, nD, all_W, all_w, tol=cfg.tol, aType=cfg.aType)
            say("%d symmetrically-unique configuration(s) for %d vacancies/dopants; time: %.2f sec." % (len(confs), nD, time.time() - initT1))
            for conf, deg in zip(confs, degens):
                # the unused interstices (X) are kept in the structure passed on to the next step, and removed from the output
                full, new_list = apply_dopants(atoms, conf, cfg.dType, eList=dexcl, dexcl_type=cfg.dexcl_type, atoms_list=atoms_list, mods=mods)
                new = full
                if cfg.aType and 'X' in cfg.aType and 'X' in full.get_chemical_symbols():
                    new = full.copy()
                    del new[[atom.index for atom in new if atom.symbol == 'X']]
                configs.append(ConfigResult(step=step, input_index=find, nD=nD, sites=list(conf) + [e for e in dexcl if e not in conf],
                                            degeneracy=int(deg), atoms=new, atoms_list=new_list, next_atoms=full))
        return configs

    def structure_keys(self, configs, root):
        """Keys of the configurations derived from the same Step 1 structure; the same for equivalent structures (see
        enumeration.canonical_keys). Only keys from the same call can be compared."""
        parent = self._parents[root]
        all_W, all_w = self._ref_symmetry[root]
        changes = [changed_sites(parent.get_chemical_symbols(), c.atoms) for c in configs]
        keys = canonical_keys(changes, parent.get_scaled_positions(wrap=True), all_W, all_w, tol=self.cfg.tol)
        return [(root, c.atoms.get_chemical_formula(), key) for c, key in zip(configs, keys)]

    def merge_equivalent(self, configs, root):
        """Merge the equivalent configurations of an input structure, adding their degeneracies.

        In the later steps of a sequential doping, the configurations are enumerated with the symmetry of the doped
        structure; two of them can still be equivalent through an operation of the parent (Step 1) structure that moves
        the dopants added earlier.
        """
        seen = {}
        out = []
        for c, key in zip(configs, self.structure_keys(configs, root)):
            if key in seen:
                seen[key].degeneracy += c.degeneracy
                continue
            seen[key] = c
            out.append(c)
        if len(out) < len(configs):
            say("%d configuration(s) equivalent to another one by the symmetry of the Step 1 structure were merged." % (len(configs) - len(out)))
        return out

    def write_configs(self, configs, outdir):
        """Write the structures of the configurations and set their names and files."""
        cfg = self.cfg
        ext = cfg.otype.split("-")[-1]  # Outfile extension.
        idir = cfg.idir if os.path.isabs(cfg.idir) else os.path.join(self.base, cfg.idir)
        copyinp = cfg.copyinp
        if copyinp and (not os.path.exists(idir) or len(os.listdir(idir)) == 0):
            say("Copying input files is requested, but the %s directory is empty. Switching off -copy." % cfg.idir)
            copyinp = False
        for cnt, conf in enumerate(configs):
            conf.name = '%s_%d' % (seed, cnt)
            if cfg.sepfol:
                d = os.path.join(outdir, str(cnt))
                os.makedirs(d, exist_ok=True)
                if cfg.outf:
                    fname = os.path.join(d, cfg.outf)
                else:
                    fname = os.path.join(d, "%s_%d.%s" % (seed, cnt, ext))
                write_struct(fname, conf.atoms, cfg.otype)

                if copyinp:  # Copy input files from INPUTS dir.
                    for fn in os.listdir(idir):
                        if os.path.isfile(os.path.join(idir, fn)):
                            shutil.copy(os.path.join(idir, fn), d)
                    if cfg.makepotcar and cfg.otype == "vasp":
                        if not os.path.exists(os.path.join(d, 'POSCAR')):
                            write_struct(os.path.join(d, 'POSCAR'), conf.atoms, "vasp")
                        make_potcar(xc=cfg.potcarxc, wDir=d)
            else:  # if all out files in the same folder.
                name = cfg.outf.split('_')[0] if cfg.outf else seed
                fname = os.path.join(outdir, "%s_%d.%s" % (name, cnt, ext))
                write_struct(fname, conf.atoms, cfg.otype)
            conf.file = os.path.relpath(fname, self.base)

    def compute(self, configs, outdir, stepdir):
        """Compute the energies and space groups of the configurations, and write the results as .res files."""
        cfg = self.cfg
        ntasks = cfg.ntasks
        if len(configs) < ntasks:
            say('Fewer configs to calculate than the set no of parallel tasks (--ntasks), reducing no of tasks to 1... ')
            ntasks = 1
        for i, E, atoms, P, err in run_calculations(self.engine, configs, outdir, ntasks=ntasks):
            conf = configs[i]
            if err is not None:
                if cfg.noskip:
                    raise RuntimeError('Problem with the calculation of %s (%s), stopping as requested...' % (conf.file, err))
                warn('Problem with the calculation of %s (%s)! \nSkipping the config and assign it with 0 energy and P1 symmetry '
                     'as placeholder...' % (conf.file, err))
                conf.failed = True
                conf.energy = 0.0
                conf.spacegroup = 'P1 (1)'
                conf.final_atoms = atoms
                continue
            conf.energy = E
            conf.pressure = P
            conf.final_atoms = atoms
            conf.spacegroup = spacegroup_string(atoms, cfg.tol)
            stoich = atoms.get_chemical_formula('hill', empirical=True)
            fn1 = os.path.join(stepdir, '%s-%s-%s.res' % (stoich, 'config_enum', hashid()))
            write_res_result(fn1, atoms, E, pressure=P, spacegroup=conf.spacegroup.split(' (')[0], name=stoich)
            conf.result_file = os.path.relpath(fn1, self.base)
            debug('%s: Ep= %.5f P= %.2f GPa SG= %s' % (conf.name, E, P, conf.spacegroup))

    def thermodynamics(self, configs):
        """Reduced energies and probabilities of the configurations (separately for each number of dopants)."""
        cfg = self.cfg
        for nD in sorted({c.nD for c in configs}):
            group = [c for c in configs if c.nD == nD]
            # probability and configurational free energy.
            P, Em_bar = boltz_dist([c.energy for c in group], T=cfg.temp, omega=[c.degeneracy for c in group], verb=cfg.verb)
            for c, p, e in zip(group, P, Em_bar):
                c.prob = float(p)
                c.E_reduced = float(e)

    def table(self, configs, inpf):
        etype = 'E_DFT' if self.engine.is_dft else 'E_%s' % self.engine.label
        str1 = "%s\n" % inpf
        str1 += "Config-ID   Degens      %-5s       E_reduced  Prob.    Space Group\n" % etype + 80 * '-' + "\n"
        for i, c in enumerate(configs):
            str1 += " config_%-3d   %-4d  %12.5f  %12.5f  %.3f    %s \n" % (i, c.degeneracy, c.energy, c.E_reduced, c.prob, c.spacegroup)
        return str1

    def sort_configs(self, configs):
        """Sort the configurations of a step for picking (-stype)."""
        seqtype = self.cfg.seqtype
        if seqtype == 'S':  # symmetry-based sorting (highest symmetry first; lowest energy first for the same symmetry)
            return sorted(configs, key=lambda c: (-spacegroup_number(c.spacegroup), c.E_reduced))
        elif seqtype == 'E':  # energy-based sorting
            return sorted(configs, key=lambda c: c.E_reduced)
        elif seqtype == 'D':  # degeneracy-based sorting
            return sorted(configs, key=lambda c: -c.degeneracy)
        elif seqtype == 'R':  # randomise
            ind = list(range(len(configs)))
            np.random.shuffle(ind)
            return [configs[int(i)] for i in ind]
        raise ValueError("Unknown sorting type (-stype): %s" % seqtype)

    def pick(self, configs, stepdir, step):
        """Pick the lucky configs for the next round; returns (picked, text)."""
        cfg = self.cfg
        x = self.sort_configs(configs)
        for root in sorted({c.root for c in x}):  # keys to recognise equivalent configs made from different input structures
            group = [c for c in x if c.root == root]
            for c, key in zip(group, self.structure_keys(group, root)):
                c.key = key
        str1 = "\nAll sorted configs for Step %d\n" % (step)
        str1 += '        Config-ID               Energy     SG \n%s\n' % (60 * '-')
        for c in x:
            str1 += '%s %s %s\n' % (c.file, c.E_reduced, spacegroup_number(c.spacegroup))
        str1 += '\n'

        # Check for duplicate configs generated from different input files at the current step.
        picked = []
        pKeys = []
        pLines = []  # Lines of previous configs
        pdir = os.path.join(stepdir, 'Picked')
        shutil.rmtree(pdir, ignore_errors=True)
        os.makedirs(pdir)
        open(os.path.join(pdir, 'Structures-are-not-optimised'), 'w').close()
        for cf in x:  # loop over ordered configs.
            if len(picked) == cfg.seqpick:
                break
            if cf.failed and not all(c.failed for c in x):
                continue
            with open(os.path.join(self.base, cf.file)) as f:
                ln1 = f.readlines()
            flag = 0
            for j, of in enumerate(pLines):  # Checking the same (or symmetrically-equivalent) entries here.
                if ln1 == of or (cf.key is not None and cf.key == pKeys[j]):
                    str1 += '%s matches to %s, skipping...\n' % (cf.file, picked[j].file)
                    flag = 1
                    break
            if flag == 1:
                continue
            cf.picked = True
            picked.append(cf)
            pKeys.append(cf.key)
            pLines.append(ln1)
            parts = cf.file.split(os.sep)
            shutil.copy(os.path.join(self.base, cf.file), os.path.join(pdir, "-".join(parts[1:])))

        str1 += '\nPicked configs for Step %d: \n%s\n' % (step, '  '.join(c.file for c in picked))
        return picked, str1

    def run_input(self, inp, step, find, nDops, stepdir):
        """All of the work for one input structure of a step; returns (configs, text of the table) or None (dry run)."""
        cfg = self.cfg
        say('Working on %s...' % inp.path)
        outdir = os.path.join(stepdir, 'input_%d' % (find))
        if os.path.exists(outdir):
            if cfg.overwrite:
                say("%s folder exists. Overwriting as requested." % os.path.relpath(outdir, self.base))
            elif cfg.restart:
                say("%s folder exists. Restarting as requested." % os.path.relpath(outdir, self.base))
            else:
                raise RuntimeError("%s folder exists. Not overwriting, use -ow to overwrite." % os.path.relpath(outdir, self.base))
        if not cfg.dry:
            os.makedirs(outdir, exist_ok=True)

        atoms, rpt_no = self.prepare_structure(inp, step, find)
        root = find if inp.root is None else inp.root
        all_W, all_w, wyckoffs = self.get_symmetry(atoms, step, root, inpf=inp.path, rpt_no=rpt_no)
        if PARENT_ARRAY not in atoms.arrays:  # a new (Step 1) structure: the parent of all the configurations made from it
            root = find
            set_parent_sites(atoms)
            self._parents[root] = atoms.copy()
            self._ref_symmetry[root] = (all_W, all_w)
        if cfg.dry:
            say("Dry run: No inputs were created.")
            return None

        atoms_list, eList = self.select(atoms, wyckoffs, inp, step)

        say("\nNow creating the modified cells, applying the desired type and number of dopant/vacancy on the selected atoms.")
        configs = self.generate(atoms, atoms_list, eList, nDops, all_W, all_w, step, find)
        if len(configs) == 0:
            raise RuntimeError("No configurations were generated. Terminating...")
        for c in configs:
            c.root = root
        if step > 1 and cfg.seqsymm == 'current':
            configs = self.merge_equivalent(configs, root)
        if cfg.ASE_symm:
            configs = ase_symmetry_filter(configs)
        say("\nTotal of %d symmetrically-unique configurations were generated." % len(configs))

        self.write_configs(configs, outdir)
        # Now do the actual calculation (SP or geom opt).
        say("Energies with: %s" % self.engine.describe())
        self.compute(configs, outdir, stepdir)
        self.thermodynamics(configs)
        return configs, self.table(configs, inp.path)

    def run(self, inpf):
        """Run the workflow on the input file(s); returns WorkflowResults.

        Args:
            inpf: a file name or a list of file names (or InputStructure objects).
        """
        cfg = self.cfg
        initT0 = time.time()
        if isinstance(inpf, (str, InputStructure)):
            inpf = [inpf]
        inputs = [i if isinstance(i, InputStructure) else InputStructure(path=i) for i in inpf]

        nDops = sorted(int(n) for n in cfg.noDopants)
        if cfg.seq > 0:
            if len(nDops) > 1:
                raise ValueError('Multiple no of dopant entries (-nd) is not supported with -seq option...')
            tDop = nDops[0]
            plan = [[min(cfg.seq, tDop - cDop)] for cDop in range(0, tDop, cfg.seq)]
        else:
            plan = [nDops]
        if not cfg.dry:
            self._check_existing()

        # Start the sequential doping loop here
        for step, step_nDops in enumerate(plan, start=1):
            initT = time.time()
            stepdir = self.step_dir(step)
            say('Starting Step %d...' % (step))
            if not cfg.dry:
                os.makedirs(stepdir, exist_ok=True)
                for fn in glob.glob(os.path.join(stepdir, '*-config_enum-??????.res')):  # results of an earlier run (-rs)
                    os.remove(fn)
            str1 = cfg.command + '\n\n' if step == 1 and cfg.command else ''
            sres = StepResult(step=step, nDops=list(step_nDops), inputs=[i.path for i in inputs])

            # Start input file loop
            for find, inp in enumerate(inputs):
                out = self.run_input(inp, step, find, step_nDops, stepdir)
                if out is None:  # dry run
                    return self.results
                configs, table = out
                str1 += table
                sres.configs.extend(configs)
            say('End of the file loop for Step %d\n' % step)

            # Now pick the lucky configs for the next round...
            if cfg.verb:
                say("All configs for Step %d" % step)
                for c in sres.configs:
                    say(c.file, c.E_reduced, spacegroup_number(c.spacegroup))
            sres.picked, text = self.pick(sres.configs, stepdir, step)
            str1 += text
            if self.engine.is_dft or self.engine.settings.optimise:
                open(os.path.join(stepdir, 'Structures-are-optimised'), 'w').close()

            sres.summary = str1
            with open(os.path.join(stepdir, 'degeneracies'), 'w') as degen_outf:
                degen_outf.write(str1 + '\n\n')
            say(str1)
            self.results.steps.append(sres)
            self.write_results()

            say('Done with Step %d...' % step)
            say("Elapsed time for Step %d: %.2f sec.\n" % (step, time.time() - initT))
            # the picked configs are the input structures of the next step
            inputs = [InputStructure(path=c.file, atoms=c.next_atoms, atoms_list=c.atoms_list, root=c.root) for c in sres.picked]
            # End of sequential loop

        say("Total elapsed time: %.2f sec." % (time.time() - initT0))
        return self.results

    def write_results(self, fname='results.json'):
        with open(os.path.join(self.base, fname), 'w') as f:
            json.dump(self.results.as_dict(), f, indent=1)
