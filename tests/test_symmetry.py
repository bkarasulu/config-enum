import numpy as np
import pytest

from config_enum.symmetry import (find_equisites, find_interstices, find_primitive_cell, get_equivalent_sites, get_symmetry,
                                  identity_symmetry, invariant_symops, read_wyckoff_data, site_permutations, spacegroup_number,
                                  spacegroup_string)

WYCKOFF = """# test data in the format of wyckoff.dat (only two space groups)
221 1 P m -3 m Oh^1 cubic
1 a m-3m
(0,0,0)
1 b m-3m
(1/2,1/2,1/2)
3 c 4/mm.m
(0,1/2,1/2) (1/2,0,1/2) (1/2,1/2,0)
3 d 4/mm.m
(1/2,0,0) (0,1/2,0) (0,0,1/2)

225 1 F m -3 m Oh^5 cubic
4 a m-3m
(0,0,0)
4 b m-3m
(1/2,1/2,1/2)
8 c -43m
(1/4,1/4,1/4) (1/4,1/4,3/4)
24 e 4m.m
(x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x)
"""


def test_spacegroup(nacl, tio2):
    assert spacegroup_string(nacl) == 'Fm-3m (225)' and spacegroup_number('Fm-3m (225)') == 225
    data = get_symmetry(tio2)
    assert data['number'] == 136 and len(data['rotations']) == 16 and len(data['translations']) == 16
    assert set(data['wyckoffs']) == {'a', 'f'}


def test_equivalent_sites(nacl, tio2):
    eSites, wyckoffs = get_equivalent_sites(nacl, get_symmetry(nacl))
    assert sorted(eSites) == ['4a', '4b']
    assert sorted(eSites['4a'] + eSites['4b']) == sorted('%s%d' % (a.symbol, a.index) for a in nacl)
    assert {wyckoffs[str(a.index)] for a in nacl if a.symbol == 'Na'} == {'4a'}
    eSites, wyckoffs = get_equivalent_sites(tio2, get_symmetry(tio2))
    assert len(eSites['2a']) == 2 and len(eSites['4f']) == 4 and wyckoffs['0'] == '2a'
    # without symmetry information every atom is its own site
    eSites, wyckoffs = get_equivalent_sites(tio2, None)
    assert len(wyckoffs) == len(tio2)


def test_primitive_cell_keeps_conventional_labels(nacl):
    prim, rpt_no = find_primitive_cell(nacl)
    assert len(prim) == 2 and rpt_no == 4
    eSites, wyckoffs = get_equivalent_sites(prim, get_symmetry(prim), rpt_no=rpt_no)
    assert sorted(eSites) == ['4a', '4b']


def test_site_permutations(nacl):
    data = get_symmetry(nacl)
    sc = nacl.get_scaled_positions(wrap=True)
    na = [a.index for a in nacl if a.symbol == 'Na']
    perms = site_permutations(sc[na], data['rotations'], data['translations'])
    assert perms.min() >= 0 and all(sorted(p) == [0, 1, 2, 3] for p in perms)
    assert len(perms) == 24  # all permutations of the four Na sites (192 operations)
    W, w = identity_symmetry()
    assert site_permutations(sc[na], W, w).tolist() == [[0, 1, 2, 3]]
    # an operation that does not map the sites onto each other is marked with -1
    perms = site_permutations(sc[na], W, np.array([[0.1, 0., 0.]]))
    assert perms.tolist() == [[-1, -1, -1, -1]]


def test_find_equisites_matches_spglib(tio2):
    data = get_symmetry(tio2)
    assert find_equisites(tio2, data['rotations'], data['translations'], tol=1e-6) == [int(i) for i in data['equivalent_atoms']]


def test_invariant_symops(nacl):
    data = get_symmetry(nacl)
    doped = nacl.copy()
    doped[[a.index for a in doped if a.symbol == 'Na'][0]].symbol = 'K'
    W, w = invariant_symops(doped, data['rotations'], data['translations'])
    assert len(W) == 48  # the point group of the K site
    own = get_symmetry(doped)
    assert len(own['rotations']) == 48


def test_read_wyckoff_data(tmp_path):
    f = tmp_path / "wyckoff.dat"
    f.write_text(WYCKOFF)
    data = read_wyckoff_data(str(f), 225)
    assert list(data) == ['225-1'] and sorted(data['225-1']) == ['24e', '4a', '4b', '8c']
    assert data['225-1']['8c'] == [[0.25, 0.25, 0.25], [0.25, 0.25, 0.75]]
    assert data['225-1']['24e'][0] == [0., 0., 0.]  # free parameters are set to 0
    assert sorted(read_wyckoff_data(str(f), 221)['221-1']) == ['1a', '1b', '3c', '3d']
    assert read_wyckoff_data(str(f), 1) == {}


def test_find_interstices(tmp_path, cu_bulk, nacl):
    f = tmp_path / "wyckoff.dat"
    f.write_text(WYCKOFF)
    cu = find_interstices(cu_bulk.copy(), tol=1e-3, wyckoff_file=str(f))
    # fcc: 4 octahedral (4b) and 8 tetrahedral (8c) interstices
    assert cu.get_chemical_symbols().count('X') == 12 and cu.get_chemical_symbols().count('Cu') == 4
    sc = cu.get_scaled_positions()[4:]
    assert sum(np.allclose(np.abs(s - np.round(s * 4) / 4), 0, atol=1e-8) for s in sc) == 12
    # rock salt: only the tetrahedral sites are free
    assert find_interstices(nacl.copy(), tol=1e-3, wyckoff_file=str(f)).get_chemical_symbols().count('X') == 8
    with pytest.raises(FileNotFoundError):
        find_interstices(cu_bulk.copy(), wyckoff_file=str(tmp_path / "missing.dat"))
