"""Atom selection, symmetry-unique configurations, doping and the final modifications (-mod)."""
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
from ase.geometry import get_distances

from .log import debug, say
from .symmetry import site_permutations
from .utils import ncr, parse_ranges, sort_dict_by_val

SITE_ARRAY = 'config_enum_site'  # temporary per-atom array: 1 for the sites still available for doping
PARENT_ARRAY = 'config_enum_parent'  # per-atom array: index of the atom in the parent (Step 1) structure


# ------------------------------------------------------------------------------------------------ atom selection
def parse_atom_list(text, noAtoms):
    """'-nl 1,2,5,8:12' (count starts at 1; '8:' = up to the last atom) -> list of atom indices starting at 0."""
    atoms_list = []
    for y in str(text).split(','):
        y = y.strip()
        if y == '':
            continue
        try:
            if ':' in y:
                x = y.split(':')
                x0 = int(x[0])
                x1 = noAtoms if len(x[1]) == 0 else int(x[1])
                vals = list(range(x0 - 1, x1))
            else:
                x0 = x1 = int(y)
                vals = [x0 - 1]
        except ValueError:
            raise ValueError('atom list(count starts at 1); e.g. -nl 1,2,5,8:12')
        if x0 <= 0 or x1 > noAtoms:
            raise ValueError("Wrong atom index (atom count starts at 1)")
        atoms_list.extend(vals)
    return atoms_list


def parse_exclude_list(excludelist):
    """'-excl 1:4 24' -> [0, 1, 2, 3, 23] (atom numbering starts at 1 on input, at 0 on output)."""
    if not excludelist:
        return []
    return sorted(x - 1 for x in parse_ranges(excludelist))


def select_atoms(atoms, aType=None, wycklist=None, wyckoffs=None, atom_list=None):
    """Indices (starting at 0) of the atoms/crystal sites to be considered in the enumeration.

    Args:
        aType: atom types (-at); 'Va' selects the X (interstice) atoms.
        wycklist: Wyckoff labels (-wl), as in the wyckoffs dict {'atom index': 'label'} of get_equivalent_sites.
            If both aType and wycklist are given, an atom must satisfy both.
        atom_list: text of the -nl option (count starts at 1); used if neither aType nor wycklist is given.
    """
    wyckoffs = wyckoffs or {}
    atoms_list = []
    if aType or wycklist:  # Determine the atomic indices of the atoms with the desired type.
        chemsym = atoms.get_chemical_symbols()
        aType = ['X' if a == 'Va' else a for a in aType] if aType else aType
        for i in range(len(chemsym)):
            at = chemsym[i]
            wt = wyckoffs.get(str(i), '')
            if aType and wycklist:
                if at in aType and wt in wycklist:
                    atoms_list.append(i)
            elif aType:
                if at in aType:
                    atoms_list.append(i)
            elif wycklist:
                if wt in wycklist:
                    atoms_list.append(i)
    elif atom_list:
        atoms_list = parse_atom_list(atom_list, len(atoms))
    return atoms_list


# ---------------------------------------------------------------------------------- symmetry-unique configurations
def get_combinations(atoms, atoms_list, nD, aType=None):
    """All combinations of nD of the selected atoms.

    If more than one atom type is selected (-at Li Nb), only the configs where all of the target types are being doped
    are included (for instance, if you choose Li and Nb sites to dope, then the configs where both Li and Nb sites
    are doped).
    """
    c = [x for x in combinations(atoms_list, r=nD)]  # gives the selected atom ids
    if aType and len(aType) > 1:
        syms = atoms.get_chemical_symbols()
        c = [x for x in c if len(np.unique([syms[i] for i in x])) == len(aType)]
    return c


def find_unique_configs(c, scpos, all_W, all_w, tol=1e-3):
    """Reduce a list of configurations (combinations of atom indices) to the symmetrically-inequivalent ones.

    Two configurations are equivalent if a symmetry operation maps the sites of one onto the sites of the other
    (in any order). The first configuration found of each set is kept.

    Args:
        c: list of tuples of atom indices.
        scpos: scaled positions of all atoms.
        all_W, all_w: rotations and translations of the symmetry operations.
        tol: tolerance on the squared distance of scaled positions.

    Returns:
        (unique, degens): indices (in c) of the irreducible configurations and the degeneracy (number of equivalent
        configurations in c) of each.
    """
    if len(c) == 0:
        return [], []
    sites = sorted({int(i) for conf in c for i in conf})
    local = {at: i for i, at in enumerate(sites)}
    perms = site_permutations(np.asarray(scpos)[sites], all_W, all_w, tol=tol)

    seen = {}  # configuration (sorted tuple of local site ids) -> index of its irreducible configuration
    unique = []
    degens = []
    for cind, conf in enumerate(c):  # iterate over configurations.
        key = tuple(sorted(local[int(at)] for at in conf))
        i = seen.get(key)
        if i is None:  # a new config: store all of its symmetry-equivalent images
            i = len(unique)
            unique.append(cind)
            degens.append(0)
            debug('New config: ', cind)
            seen[key] = i
            if len(key):
                images = perms[:, list(key)]
                images = images[np.all(images >= 0, axis=1)]
                for img in np.unique(np.sort(images, axis=1), axis=0):
                    seen.setdefault(tuple(int(x) for x in img), i)
        degens[i] += 1
    return unique, degens


# ------------------------------------------------------------------------------------------------- modifications
@dataclass
class Modification:
    """A final modification of the configurations (-mod), e.g. 'replace all within 2 of Al with Va'."""
    action: str = ''          # 'replace' or 'remove'
    mod_repl: list = field(default_factory=list)  # atom types that can be replaced/removed
    mod_type: str = ''        # 'within' or 'nearest'
    mod_value: float = 0.     # distance (A) for 'within'; number of neighbours for 'nearest'
    mod_target: str = ''      # atom type used as reference, or 'dopant'
    mod_replwith: str = None  # new atom type for 'replace' ('Va': delete)


def parse_modify(text, atomSyms=()):
    """Parse the --modify keyword: '<replace|remove> <atom type|all> <within|nearest> <value> of <target> [with <new>]'.

    Several statements can be given, separated by ';'; they are applied in the given order.
    Returns a list of Modification objects.
    """
    mods = []
    if not text:
        return mods
    for statement in str(text).replace("'", "").replace('"', '').split(';'):
        mod_keys = statement.split()
        if len(mod_keys) == 0:
            continue
        flag = 0
        if len(mod_keys) < 6:
            flag = 1
        elif mod_keys[0] not in ["remove", "replace"]:
            flag = 1
        elif mod_keys[2] not in ["within", "nearest"]:
            flag = 1
        elif mod_keys[4] not in ["of"]:
            flag = 1
        elif mod_keys[0] == "replace" and (len(mod_keys) != 8 or mod_keys[6] != "with"):
            flag = 1
        if flag:
            raise ValueError("wrong usage of -mod keyword: '%s'. Example: -mod 'replace all within 2 of Al with Va'" % statement.strip())

        mod = Modification(action=mod_keys[0])
        mod.mod_repl = list(atomSyms) if mod_keys[1] == "all" else [mod_keys[1]]
        mod.mod_type = mod_keys[2]  # within or nearest
        mod.mod_value = float(mod_keys[3])
        mod.mod_target = mod_keys[5]
        if mod.action == "replace":
            mod.mod_replwith = mod_keys[-1]
        mods.append(mod)
    return mods


def apply_modification(atoms, mod, dopant_positions=()):
    """Apply a final modification (-mod) to a configuration (in place).

    Args:
        mod: Modification
        dopant_positions: Cartesian positions of the dopants/vacancies of this configuration (reference points if the
            target is 'dopant').
    Returns the atoms.
    """
    if mod.mod_type not in ['within', 'nearest']:
        return atoms
    search_list = [atom.index for atom in atoms if atom.symbol in mod.mod_repl]  # atoms that can be modified
    # this is the list of positions used as reference (i.e. mod_target).
    if mod.mod_target.lower() in ['do', 'dopant', 'dopants']:
        ref_list = [(None, np.asarray(p, dtype=float)) for p in dopant_positions]
    else:
        ref_list = [(atom.index, atom.position.copy()) for atom in atoms if atom.symbol == mod.mod_target]
    sel_list = []  # this is the list atoms that will be modified
    if len(ref_list) == 0 or len(search_list) == 0:
        return atoms

    _, D = get_distances(np.array([r[1] for r in ref_list]), atoms.positions[search_list], cell=atoms.cell, pbc=atoms.pbc)
    cn = 0
    dists0 = {}
    for n1, (ind1, _pos) in enumerate(ref_list):
        dists = {}
        for n2, ind2 in enumerate(search_list):
            if ind1 == ind2:
                continue
            if ind1 is None and D[n1, n2] < 1e-6:
                continue  # the dopant itself
            dist = D[n1, n2]
            dists0["%s-%s" % (n1, ind2)] = dist
            dists[str(ind2)] = dist
            if mod.mod_type == 'within' and dist <= float(mod.mod_value):
                if ind2 not in sel_list:
                    sel_list.append(ind2)

        if mod.mod_type == 'nearest':
            if len(dists) == 0:
                continue
            keys = sort_dict_by_val(dists)
            for key in keys[0:int(mod.mod_value)]:
                if int(key) not in sel_list:
                    sel_list.append(int(key))
        cn += 1

    # This is needed for adding the closest atom to any of the reference atoms to top up to the requested number of neighbours.
    if mod.mod_type == 'nearest' and len(sel_list) < cn * int(mod.mod_value):
        for ky in sort_dict_by_val(dists0):
            key = int(ky.split('-')[1])
            if key not in sel_list:
                sel_list.append(key)
            if len(sel_list) == cn * int(mod.mod_value):
                break

    if mod.action == "replace" and mod.mod_replwith != 'Va':
        for sel in sel_list:
            atoms[sel].symbol = mod.mod_replwith
            if SITE_ARRAY in atoms.arrays:
                atoms.arrays[SITE_ARRAY][sel] = 0
    else:  # remove, or replace with Va
        del atoms[sel_list]
    return atoms


# ---------------------------------------------------------------------------------------------------------- doping
def apply_dopants(atoms_orig, conf, dType, eList=(), dexcl_type=None, atoms_list=None, mods=(), delete_X=False):
    """Create a configuration: put the dopant (or vacancy) on the given sites of a structure.

    Args:
        atoms_orig: the original structure (not modified).
        conf: indices of the atoms to be replaced by dType ('Va': the atoms are deleted).
        eList: indices of further atoms to be replaced by dexcl_type (the excluded list with -dexcl).
        dexcl_type: dopant type for eList (default: dType).
        atoms_list: the atoms selected for the enumeration; used to work out which sites are still available in the
            new structure (for sequential doping).
        mods: list of Modification objects, applied in the given order.
        delete_X: delete the remaining X atoms (interstices added in the beginning).

    Returns:
        (atoms, atoms_list_new): the new structure and the indices (in the new structure) of the selected sites that
        were not doped or modified.
    """
    atoms = atoms_orig.copy()
    conf = [int(x) for x in conf]
    eList = [int(x) for x in eList]
    if dexcl_type in (None, ''):
        dexcl_type = dType

    site = np.zeros(len(atoms), dtype=int)
    if atoms_list is not None and len(atoms_list):
        site[list(atoms_list)] = 1
    site[conf] = 0
    site[eList] = 0
    atoms.set_array(SITE_ARRAY, site)
    dopant_positions = atoms.positions[conf + eList].copy()

    # Prepare the full coordinate set for the given config.
    to_del = []
    for at in conf:
        if dType == 'Va':
            to_del.append(at)
        else:
            atoms[at].symbol = dType
    for at in eList:
        if at in conf:
            continue
        if dexcl_type == 'Va':
            to_del.append(at)
        else:
            atoms[at].symbol = dexcl_type
    if to_del:
        del atoms[to_del]  # multiple deletion works well.

    # Add the final modifications here!!!
    for mod in mods:
        apply_modification(atoms, mod, dopant_positions=dopant_positions)

    # delete the X atoms added for determining the interstitial sites in the beginning.
    if delete_X:
        del atoms[[atom.index for atom in atoms if atom.symbol == 'X']]

    atoms_list_new = [int(i) for i in np.where(atoms.arrays[SITE_ARRAY] == 1)[0]]
    del atoms.arrays[SITE_ARRAY]
    return atoms, atoms_list_new


def enumerate_configs(atoms, atoms_list, nD, all_W, all_w, tol=1e-3, aType=None):
    """The symmetrically-unique configurations of nD dopants/vacancies on the selected atoms.

    Returns (confs, degens): list of tuples of atom indices and their degeneracies.
    """
    nc = ncr(len(atoms_list), nD)
    if nc > 2e4:  # Too many possibilities.
        say("Large number of combinations (%.1e) for %d vacancies/dopants, this can take a while..." % (nc, nD))
    c = get_combinations(atoms, atoms_list, nD, aType=aType)
    say("%d configuration(s) found for %d vacancies/dopants." % (len(c), nD))
    scpos = atoms.get_scaled_positions(wrap=True)
    unique, degens = find_unique_configs(c, scpos, all_W, all_w, tol=tol)
    return [c[i] for i in unique], degens


# ------------------------------------------------------- equivalence of whole structures (sequential doping)
def set_parent_sites(atoms):
    """Label the atoms with their index in this (parent) structure; the labels survive doping and deletions."""
    atoms.set_array(PARENT_ARRAY, np.arange(len(atoms)))
    return atoms


def changed_sites(parent_symbols, atoms):
    """How a structure differs from its parent: {parent site: new atom type ('Va' if the atom was removed)}."""
    present = dict(zip((int(i) for i in atoms.arrays[PARENT_ARRAY]), atoms.get_chemical_symbols()))
    return {s: present.get(s, 'Va') for s in range(len(parent_symbols)) if present.get(s, 'Va') != parent_symbols[s]}


def canonical_keys(changes, parent_scpos, all_W, all_w, tol=1e-3):
    """A key for each structure that is the same for symmetrically-equivalent structures.

    Each structure is given by its changes with respect to a common parent structure (see changed_sites); two
    structures are equivalent if a symmetry operation of the parent maps the changes of one onto those of the other.
    This also finds the equivalence of configurations created by adding dopants in a different order.

    Args:
        changes: list of dicts {parent site: atom type}
        parent_scpos: scaled positions of the parent structure.
        all_W, all_w: symmetry operations of the parent structure.
    """
    sites = sorted({s for ch in changes for s in ch})
    if len(sites) == 0:
        return [() for ch in changes]
    local = {s: i for i, s in enumerate(sites)}
    labels = sorted({v for ch in changes for v in ch.values()})
    code = {v: i for i, v in enumerate(labels)}
    L = len(labels)
    perms = site_permutations(np.asarray(parent_scpos)[sites], all_W, all_w, tol=tol)
    keys = []
    for ch in changes:
        if len(ch) == 0:
            keys.append(())
            continue
        idx = [local[s] for s in ch]
        lab = np.array([code[v] for v in ch.values()])
        images = perms[:, idx]
        images = images[np.all(images >= 0, axis=1)]
        arr = np.sort(images * L + lab[None, :], axis=1)
        first = arr[np.lexsort(arr.T[::-1])[0]]  # the smallest image
        keys.append(tuple((sites[int(x) // L], labels[int(x) % L]) for x in first))
    return keys
