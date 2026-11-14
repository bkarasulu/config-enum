"""Symmetry operations, equivalent (Wyckoff) sites, primitive cells and interstices (spglib and ASE)."""
import os
import warnings
from copy import deepcopy
from fractions import Fraction

import numpy as np
import spglib
from ase import Atom, Atoms

from .geometry import volume
from .log import debug, say, warn
from .utils import sorted_nicely


def _cell(atoms):
    return (np.array(atoms.cell), atoms.get_scaled_positions(), atoms.numbers)


def spacegroup_string(atoms, tol=1e-3):
    """Space group as given by spglib, e.g. 'Fm-3m (225)' ('P1 (1)' if it cannot be determined)."""
    sg = spglib.get_spacegroup(_cell(atoms), symprec=tol)
    return sg if sg else 'P1 (1)'


def spacegroup_number(sg):
    """'Fm-3m (225)' -> 225"""
    return int(str(sg).split('(')[-1].replace(')', ''))


def get_symmetry(atoms_crys, tol=1e-3):
    """Symmetry dataset of a structure from spglib.

    Returns a dict with 'rotations' (all_W), 'translations' (all_w), 'equivalent_atoms', 'wyckoffs', 'spacegroup'
    (e.g. 'Fm-3m (225)') and 'number'; or None if spglib cannot find the symmetry.
    """
    cell = _cell(atoms_crys)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = spglib.get_symmetry_dataset(cell, symprec=tol, angle_tolerance=-1.0, hall_number=0)
    if data is None:
        return None

    def get(key):
        return getattr(data, key) if hasattr(data, key) else data[key]
    return {'rotations': np.array(get('rotations')), 'translations': np.array(get('translations')),
            'equivalent_atoms': np.array(get('equivalent_atoms')), 'wyckoffs': list(get('wyckoffs')),
            'spacegroup': spacegroup_string(atoms_crys, tol), 'number': int(get('number'))}


def identity_symmetry():
    """The nosymm case: only the identity."""
    all_W = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
    all_w = np.zeros((1, 3))
    return all_W, all_w


def find_primitive_cell(atoms, tol=1e-3):
    """Primitive cell using SPGlib; returns (primitive atoms, rpt_no = volume ratio of the input and primitive cells)."""
    res = spglib.find_primitive(_cell(atoms), symprec=tol)
    if res is None:
        warn("WARNING: the primitive cell could not be found; using the input cell.")
        return atoms, 1
    lattice, scaled_positions, numbers = res
    prim = Atoms(symbols=numbers, cell=lattice, scaled_positions=scaled_positions, pbc=True)
    rpt_no = max(1, int(round(volume(atoms.cell) / volume(prim.cell))))
    return prim, rpt_no


def get_equivalent_sites(atoms, data, rpt_no=1):
    """Group the atoms into equivalent (Wyckoff) sites.

    Args:
        atoms: the structure whose atoms are labelled.
        data: symmetry dataset from get_symmetry (None: every atom is its own site).
        rpt_no: number of primitive cells in the original (conventional) cell, to label the sites of a primitive cell
            with the multiplicity of the conventional one.

    Returns:
        eSites: {'24d': ['Li0', 'Li1', ...], ...} (Wyckoff label -> atoms, index starts at 0)
        wyckoffs: {'0': '24d', ...} (atom index as str -> Wyckoff label)
    """
    # Determine the equivalent atoms and sites.
    equiSites = {}
    if data is not None:
        for i in range(len(data['equivalent_atoms'])):
            ea = int(data['equivalent_atoms'][i])
            if ea not in equiSites:
                try:
                    equiSites[ea] = [data['wyckoffs'][i], i]
                except IndexError:
                    equiSites[ea] = ['X', i]
            else:
                equiSites[ea].append(i)
    else:
        say("No equivalent positions. All selected atoms will be treated separately.")
        for i in range(len(atoms)):
            equiSites[i] = ['a', i]

    # Label each set of equivalent atoms with its multiplicity and Wyckoff letter; sets with the same label (e.g. two
    # atom types on 96h sites) are listed together.
    eSites = {}
    wyckoffs = {}
    for key in sorted(equiSites.keys()):  # loop over equi sites.
        val = equiSites[key]
        site = '%d%s' % (rpt_no * (len(val) - 1), val[0])
        for i in val[1:]:
            wyckoffs[str(i)] = site
            if i < len(atoms):
                eSites.setdefault(site, []).append('%s%d' % (atoms[i].symbol, i))
    return eSites, wyckoffs


def print_equivalent_sites(eSites, inpf=''):
    say("The following equivalent sites were extracted from %s: " % inpf)
    for eS in sorted_nicely(list(eSites.keys())):
        say("%s: %s" % (eS, sorted_nicely(eSites[eS])))
    say()


def invariant_symops(atoms, all_W, all_w, tol=1e-3):
    """The subset of the symmetry operations that leave the structure invariant (each atom is mapped onto an atom of
    the same type). Used to restrict the operations of a reference structure to those valid for a doped structure."""
    scpos = atoms.get_scaled_positions(wrap=True)
    numbers = atoms.numbers
    keep = []
    for i in range(len(all_W)):
        new = scpos @ np.array(all_W[i]).T + np.array(all_w[i]).reshape(3)
        ok = True
        for z in np.unique(numbers):
            a = new[numbers == z]
            b = scpos[numbers == z]
            d = a[:, None, :] - b[None, :, :]
            d -= np.round(d)
            if not np.all(np.min(np.sum(d * d, axis=-1), axis=1) <= tol):
                ok = False
                break
        if ok:
            keep.append(i)
    return np.array(all_W)[keep], np.array(all_w)[keep]


def site_permutations(scpos, all_W, all_w, tol=1e-3):
    """How the symmetry operations permute a set of sites.

    Args:
        scpos (n x 3): scaled positions of the sites.
        all_W, all_w: rotations and translations of the symmetry operations (in the basis of the cell).
        tol: two sites match if the squared distance of their scaled positions (periodic) is <= tol.

    Returns:
        (n_unique_ops x n) int array: perms[o, i] is the site onto which operation o maps site i (-1: no such site).
    """
    scpos = np.asarray(scpos, dtype=float).reshape(-1, 3)
    n = len(scpos)
    perms = np.full((len(all_W), n), -1, dtype=int)
    if n == 0:
        return perms[:1]
    for o in range(len(all_W)):
        new = scpos @ np.asarray(all_W[o], dtype=float).T + np.asarray(all_w[o], dtype=float).reshape(3)
        d = new[:, None, :] - scpos[None, :, :]
        d -= np.round(d)  # periodic boundary conditions (0.0 and 1.0 are equivalent)
        d2 = np.sum(d * d, axis=-1)
        j = np.argmin(d2, axis=1)
        ok = d2[np.arange(n), j] <= tol
        perms[o, ok] = j[ok]
    return np.unique(perms, axis=0)


def find_equisites(atoms, all_W, all_w, tol=1e-6):
    """Equivalent atoms from the symmetry operations (own implementation; an alternative to spglib's
    equivalent_atoms): returns for each atom the index of the first atom equivalent to it."""
    noAtoms = len(atoms)
    equiSites = np.arange(noAtoms)
    cpos = atoms.get_scaled_positions()  # use the scaled positions instead!!
    coords = np.zeros((4, 1))
    done = []  # already matched atoms.
    for at, _atom in enumerate(atoms):  # coords of the current atom.
        if at in done:
            continue
        coords[0:3, 0] = cpos[at][0:3]
        coords[3, 0] = 1.
        for i, W in enumerate(all_W):  # Find symmetry equi sites and match them to other available atoms.
            trans = np.zeros((4, 4))
            w = all_w[i]
            trans[0:3, 0:3] = W[:, :]
            trans[0:3, 3] = w[0:3]
            trans[3, :] = [0., 0., 0., 1.]
            new_coords = np.dot(trans, coords)
            # now match it to other atoms.
            for j in range(noAtoms):
                if j in done or atoms[j].symbol != atoms[at].symbol:
                    continue
                d = new_coords[0:3, 0] - cpos[j]  # coords of other atom.
                d -= np.round(d)  # assure the coordinates are comparable (0. and 1. are equivalent).
                if np.dot(d, d) < tol:
                    equiSites[j] = int(at)
                    done.append(j)

    return [int(x) for x in equiSites]


BUNDLED_WYCKOFF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'wyckoff.dat')


def wyckoff_file_path(wyckoff_file=None):
    """Location of wyckoff.dat: the argument, $CONFIG_ENUM_WYCKOFF, $HOME/wyckoff.dat (as in the original script) if it
    exists, or the file distributed with the package."""
    home = os.path.join(os.path.expanduser('~'), 'wyckoff.dat')
    return wyckoff_file or os.environ.get('CONFIG_ENUM_WYCKOFF') or (home if os.path.exists(home) else BUNDLED_WYCKOFF)


def read_wyckoff_data(fname, sgno, fixed_only=False):
    """Read the Wyckoff positions of space group `sgno` from wyckoff.dat.

    Format (as in the file taken from VESTA): a header line per space-group setting ('<number> <setting> ... <symbol>'),
    then for each Wyckoff position a line '<multiplicity> <letter> <site symmetry>' followed by lines of coordinates
    '(x,y,z) (-x,-y,z) ...'; the settings are separated by empty lines. The free parameters (x, y, z) are set to 0;
    with fixed_only, the positions with free parameters are left out.

    Returns {'<number>-<setting>': {'4a': [[0., 0., 0.], ...], ...}}
    """
    data = {}
    tmp = {}
    label = None
    wyck = None

    def store():
        if label is not None and tmp and int(label.split('-')[0]) == sgno:  # only store the relevant space group.
            data[label] = dict(tmp)

    with open(fname, 'r') as file:
        for ln in file:
            if '#' in ln[0:1]:
                continue
            x = ln.split()
            if ln.strip() == '':
                store()
                tmp = {}
                label = None
                continue
            elif ln[0] != '(' and len(x) >= 2 and x[0].isdigit() and len(x[1]) == 1 and x[1].isalpha():
                wyck = '%s%s' % (x[0], x[1])  # the site symmetry can contain spaces, e.g. '48 i m.m 2'
                tmp[wyck] = []
            elif ln[0] != '(' and len(x) > 5:
                label = '%s-%s' % (x[0], x[1])
            elif ln[0] == '(':
                for crd in x:  # TODO: Double check if the x,y,z and addition/substraction of fractions are handled correctly.
                    crd_bkp = deepcopy(crd)
                    if fixed_only and any(v in crd for v in 'xyz'):
                        tmp.pop(wyck, None)
                        continue
                    if wyck not in tmp:
                        continue
                    crd = crd.strip('(').strip(')')
                    for v in ('-x', '-y', '-z', '+x', '+y', '+z', 'x', 'y', 'z'):
                        crd = crd.replace(v, '')
                    try:
                        tmp[wyck].append([float(Fraction(a)) if a != '' else 0. for a in crd.split(',')])
                    except (ValueError, ZeroDivisionError, KeyError):
                        warn('Error! :  ', label, wyck, crd_bkp, crd)
    store()
    return data


def find_interstices(atoms, tol=1e-6, sg=None, verb=None, wyckoff_file=None):
    """Add the unoccupied Wyckoff positions of the space group as X atoms (potential interstitial sites).

    Only the positions without free parameters (e.g. the octahedral and tetrahedral sites of fcc) are added.

    TODO: only compatible with conventional cells (not with primitive cells) in the first setting of the space group
    (e.g. origin choice 1). Easiest solution: first add the interstitials and then find the primitive cell. May not
    work with the repeated supercells (even with conventional ones), start from standard cells.

    Args:
        atoms: structure (modified in place and returned).
        sg (int): space-group number (determined from the structure if not given).
        wyckoff_file: the wyckoff.dat file with the Wyckoff positions (see wyckoff_file_path).
    """
    from ase.spacegroup import Spacegroup

    out = say if verb else debug
    # get spacegroup
    if not sg:
        sg = spacegroup_number(spacegroup_string(atoms, tol))
    SG = Spacegroup(int(sg))
    out("Original space group: ", SG.no, SG.symbol)

    scpos = atoms.get_scaled_positions()
    out('Original scaled/fractional positions: \n', scpos)

    fname = wyckoff_file_path(wyckoff_file)
    if not os.path.exists(fname):
        raise FileNotFoundError("The Wyckoff-position data file (%s) was not found; it is needed for adding the "
                                "interstices (-ai). Give it with --wyckoff_file or $CONFIG_ENUM_WYCKOFF." % fname)
    out('\nReading %s...' % fname)
    # Only the Wyckoff positions without free parameters are used: a position such as (x,x,x) is a line of possible
    # sites, and any single value of x (the original script used 0) would be an arbitrary choice.
    data = read_wyckoff_data(fname, SG.no, fixed_only=True)

    data2 = {}
    for key in list(data.keys()):
        data2[key] = {}
        for w in list(data[key].keys()):
            if len(data[key][w]) == 0:
                continue
            newcrd = data[key][w][0]
            equiSites = SG.equivalent_sites(newcrd, onduplicates='replace')[0]  # Equi site of the first site should cover the other sites as well.
            data2[key][w] = []
            for crd in equiSites:
                flag = 0
                for sc in scpos:
                    d = np.asarray(sc) - np.asarray(crd)
                    d -= np.round(d)
                    if np.allclose(d, 0., atol=1e-5):
                        flag = 1
                        break
                if not flag:
                    data2[key][w].append(crd)

    a, b, c = atoms.get_cell()
    Va_sites = []
    for i, key in enumerate(sorted_nicely(list(data2.keys()))):
        if i == 1:
            break  # Take only the first setting
        out(key)
        sites = []
        for w in sorted_nicely(list(data2[key].keys())):
            if len(data2[key][w]) != 0:
                out(w, [" ".join(['%.3f' % y for y in x]) for x in data2[key][w]])
                sites.extend(data2[key][w])

        ulist = np.unique(np.array(sites), axis=0) if len(sites) else np.array([])

        # Sometimes np.unique fails to detect equal arrays/sites, so this is a fail-safe approach to delete repeating sites.
        for crd in ulist:
            flag = 0
            for u in Va_sites:
                if np.allclose(crd, u):
                    flag = 1
                    break
            if flag:
                continue
            r = crd[0] * a + crd[1] * b + crd[2] * c  # fractional to Cartesian conversion, works!!
            atoms.append(Atom('X', r))
            Va_sites.append(crd)

    out('\n%d unique vacancy sites detected.' % (len(Va_sites)))
    for x in Va_sites:
        out(x)
    say("\nNew space group after adding all interstices: ", spacegroup_string(atoms, tol))
    out('New scaled/fractional positions: \n', atoms.get_scaled_positions())
    return atoms
