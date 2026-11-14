from itertools import combinations

import numpy as np
import pytest
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from legacy_reference import legacy_unique_configs

from config_enum.enumeration import (apply_dopants, canonical_keys, changed_sites, enumerate_configs, find_unique_configs,
                                     get_combinations, parse_atom_list, parse_exclude_list, parse_modify, select_atoms,
                                     set_parent_sites)
from config_enum.symmetry import get_equivalent_sites, get_symmetry, identity_symmetry


def _sym(atoms):
    d = get_symmetry(atoms)
    return d['rotations'], d['translations']


def test_selection(nacl):
    _, wyckoffs = get_equivalent_sites(nacl, get_symmetry(nacl))
    na = [a.index for a in nacl if a.symbol == 'Na']
    assert select_atoms(nacl, aType=['Na']) == na
    assert select_atoms(nacl, wycklist=['4a'], wyckoffs=wyckoffs) == na
    assert select_atoms(nacl, aType=['Cl'], wycklist=['4a'], wyckoffs=wyckoffs) == []
    assert select_atoms(nacl, atom_list='1,2,5:7') == [0, 1, 4, 5, 6]
    assert parse_atom_list('7:', 8) == [6, 7]
    for bad in ('0,1', '9', 'a'):
        with pytest.raises(ValueError):
            parse_atom_list(bad, 8)
    assert parse_exclude_list(['1:4', '24']) == [0, 1, 2, 3, 23]


def test_combinations_with_two_atom_types(nacl):
    c = get_combinations(nacl, list(range(8)), 2, aType=['Na', 'Cl'])
    assert len(c) == 16  # only the configs where both types are doped
    assert len(get_combinations(nacl, list(range(8)), 2, aType=['Na'])) == 28


@pytest.mark.parametrize("rep,species,nD,expected", [((1, 1, 1), 'Na', 2, 1), ((2, 1, 1), 'Na', 2, 4), ((2, 1, 1), 'Na', 3, 4)])
def test_unique_configs_nacl(nacl, rep, species, nD, expected):
    atoms = nacl * rep
    sel = [a.index for a in atoms if a.symbol == species]
    W, w = _sym(atoms)
    c = list(combinations(sel, nD))
    unique, degens = find_unique_configs(c, atoms.get_scaled_positions(wrap=True), W, w)
    assert len(unique) == expected
    assert sum(degens) == len(c)
    # the same result as the algorithm of the original script (with the positions of the selected atoms)
    assert (unique, degens) == legacy_unique_configs(c, atoms.get_scaled_positions(wrap=True), W, w, shift=0)


def test_unique_configs_are_inequivalent_and_complete(tio2):
    """Check with ASE's structure comparison: no two unique configs are equivalent, every other config has a partner."""
    atoms = tio2 * (2, 1, 2)
    sel = [a.index for a in atoms if a.symbol == 'O']
    W, w = _sym(atoms)
    confs, degens = enumerate_configs(atoms, sel, 2, W, w)
    assert sum(degens) == 120 and len(confs) == 10
    comp = SymmetryEquivalenceCheck()
    structs = [apply_dopants(atoms, c, 'F')[0] for c in confs]
    for i in range(len(structs)):
        for j in range(i):
            assert not comp.compare(structs[i], structs[j])
    other, _ = apply_dopants(atoms, (sel[3], sel[9]), 'F')
    assert any(comp.compare(other, s) for s in structs)


def test_degeneracy_is_orbit_size(cu_bulk):
    atoms = cu_bulk * (2, 2, 1)
    W, w = _sym(atoms)
    confs, degens = enumerate_configs(atoms, list(range(16)), 1, W, w)
    assert len(confs) == 1 and degens == [16]
    confs, degens = enumerate_configs(atoms, list(range(16)), 0, W, w)
    assert confs == [()] and degens == [1]
    W, w = identity_symmetry()
    confs, degens = enumerate_configs(atoms, list(range(16)), 2, W, w)
    assert len(confs) == 120 and set(degens) == {1}


def test_apply_dopants(nacl):
    na = [a.index for a in nacl if a.symbol == 'Na']
    new, left = apply_dopants(nacl, na[:1], 'K', atoms_list=na)
    assert new.get_chemical_formula() == 'Cl4KNa3' and left == na[1:]
    assert nacl.get_chemical_formula() == 'Cl4Na4'  # the original is not changed
    new, left = apply_dopants(nacl, na[:2], 'Va', atoms_list=na)
    assert new.get_chemical_formula() == 'Cl4Na2' and len(left) == 2
    assert all(new[i].symbol == 'Na' for i in left)
    np.testing.assert_allclose(new.positions[left], nacl.positions[na[2:]])
    # excluded list with its own dopant type
    new, left = apply_dopants(nacl, na[:1], 'Va', eList=[na[3]], dexcl_type='K', atoms_list=na)
    assert new.get_chemical_formula() == 'Cl4KNa2' and len(left) == 2
    new, left = apply_dopants(nacl, na[:1], 'K', eList=[na[3]], dexcl_type='Va', atoms_list=na)
    assert new.get_chemical_formula() == 'Cl4KNa2'
    new, left = apply_dopants(nacl, na[:1], 'K', eList=[na[3]], atoms_list=na)  # default: the same type
    assert new.get_chemical_formula() == 'Cl4K2Na2'


def test_delete_interstices(cu_bulk):
    from ase import Atom
    atoms = cu_bulk.copy()
    atoms.append(Atom('X', (1.8, 1.8, 1.8)))
    atoms.append(Atom('X', (0.9, 0.9, 0.9)))
    new, left = apply_dopants(atoms, [4], 'H', atoms_list=[4, 5], delete_X=True)
    assert sorted(new.get_chemical_symbols()) == ['Cu'] * 4 + ['H'] and left == []


def test_parse_modify():
    mods = parse_modify("'replace all within 2 of Al with Va'", atomSyms=['Al', 'Li', 'O'])
    assert len(mods) == 1
    m = mods[0]
    assert (m.action, m.mod_repl, m.mod_type, m.mod_value, m.mod_target, m.mod_replwith) == ('replace', ['Al', 'Li', 'O'], 'within', 2.0, 'Al', 'Va')
    mods = parse_modify("remove Li nearest 1 of dopant; replace O within 2.5 of Li with F")
    assert [m.action for m in mods] == ['remove', 'replace'] and mods[0].mod_replwith is None
    assert parse_modify(None) == []
    for bad in ("delete Li within 2 of Al", "remove Li around 2 of Al", "replace Li within 2 of Al", "remove Li within 2 Al"):
        with pytest.raises(ValueError):
            parse_modify(bad)


def test_modifications(nacl):
    atoms = nacl * (2, 2, 2)
    na = [a.index for a in atoms if a.symbol == 'Na']
    # Mg on a Na site, and the nearest Na removed for charge balance
    new, left = apply_dopants(atoms, na[:1], 'Mg', atoms_list=na, mods=parse_modify("remove Na nearest 1 of Mg"))
    assert new.get_chemical_formula() == 'Cl32MgNa30' and len(left) == 30
    mg = [a.index for a in new if a.symbol == 'Mg'][0]
    d = np.sort(new.get_distances(mg, [a.index for a in new if a.symbol == 'Na'], mic=True))
    assert len(d[d < 5.64 / np.sqrt(2) + 0.01]) == 11  # one of the twelve nearest Na neighbours is gone
    # the same with the dopant as the reference
    new2, _ = apply_dopants(atoms, na[:1], 'Mg', atoms_list=na, mods=parse_modify("remove Na nearest 1 of dopant"))
    assert new2.get_chemical_formula() == 'Cl32MgNa30'
    # a vacancy as the reference: the six Cl neighbours within 3 A become F
    new, _ = apply_dopants(atoms, na[:1], 'Va', atoms_list=na, mods=parse_modify("replace Cl within 3 of dopant with F"))
    assert new.get_chemical_formula() == 'Cl26F6Na31'
    # replace with Va deletes
    new, _ = apply_dopants(atoms, na[:1], 'Mg', atoms_list=na, mods=parse_modify("replace Cl within 3 of Mg with Va"))
    assert new.get_chemical_formula() == 'Cl26MgNa31'
    # 'all' atom types and several statements
    mods = parse_modify("replace all within 3 of Mg with Br; remove Na nearest 2 of Mg", atomSyms=['Cl', 'Na'])
    new, left = apply_dopants(atoms, na[:1], 'Mg', atoms_list=na, mods=mods)
    assert new.get_chemical_formula() == 'Br6Cl26MgNa29' and len(left) == 29


def test_canonical_keys_find_order_independent_equivalence(cu_bulk):
    """Two dopants added one after the other: A then B is the same structure as B then A, and symmetry images too."""
    parent = set_parent_sites(cu_bulk * (2, 2, 1))
    W, w = _sym(parent)
    sym = parent.get_chemical_symbols()
    sc = parent.get_scaled_positions(wrap=True)

    def make(*sites):
        atoms = parent
        for s in sites:
            idx = [i for i, p in enumerate(atoms.arrays['config_enum_parent']) if p == s]
            atoms, _ = apply_dopants(atoms, idx, 'Va')
        return atoms
    structs = [make(0, 5), make(5, 0), make(1, 4), make(0, 1), make(0, 2), make(0)]
    changes = [changed_sites(sym, s) for s in structs]
    assert changes[0] == {0: 'Va', 5: 'Va'} and changes[5] == {0: 'Va'}
    keys = canonical_keys(changes, sc, W, w)
    comp = SymmetryEquivalenceCheck()
    for i in range(5):
        for j in range(i):
            assert (keys[i] == keys[j]) == bool(comp.compare(structs[i], structs[j])), (i, j)
    assert keys[0] == keys[1] and keys[5] != keys[0]
    assert canonical_keys([{}, {}], sc, W, w) == [(), ()]
