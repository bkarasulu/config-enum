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
    assert sorted(read_wyckoff_data(str(f), 225, fixed_only=True)['225-1']) == ['4a', '4b', '8c']
    assert sorted(read_wyckoff_data(str(f), 221)['221-1']) == ['1a', '1b', '3c', '3d']
    assert read_wyckoff_data(str(f), 1) == {}


def test_bundled_wyckoff_data():
    """The wyckoff.dat distributed with the package has all 230 space groups."""
    from config_enum.symmetry import BUNDLED_WYCKOFF, wyckoff_file_path
    assert wyckoff_file_path("my.dat") == "my.dat"
    for no in range(1, 231):
        data = read_wyckoff_data(BUNDLED_WYCKOFF, no)
        assert data and all(len(v) > 0 for d in data.values() for v in d.values()), no
    d = read_wyckoff_data(BUNDLED_WYCKOFF, 225)['225-1']
    assert sorted(d) == ['192l', '24d', '24e', '32f', '48g', '48h', '48i', '4a', '4b', '8c', '96j', '96k']
    assert d['48i'][0] == [0.5, 0., 0.] and d['24d'][0] == [0., 0.25, 0.25]  # site symmetries with a space ('m.m 2')
    assert sorted(read_wyckoff_data(BUNDLED_WYCKOFF, 225, fixed_only=True)['225-1']) == ['24d', '4a', '4b', '8c']
    assert len(read_wyckoff_data(BUNDLED_WYCKOFF, 227)) == 2  # two origin choices


def _sites(atoms):
    eSites, _ = get_equivalent_sites(atoms, get_symmetry(atoms))
    return {k: len(v) for k, v in eSites.items() if v[0].startswith('X')}


def test_find_interstices(tmp_path, cu_bulk, nacl, tio2, monkeypatch):
    from ase.build import bulk
    monkeypatch.delenv("CONFIG_ENUM_WYCKOFF", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))  # no wyckoff.dat there: the file of the package is used
    cu = find_interstices(cu_bulk.copy(), tol=1e-3)
    assert spacegroup_string(cu) == 'Fm-3m (225)'
    assert _sites(cu) == {'4b': 4, '8c': 8, '24d': 24}  # octahedral, tetrahedral and (0,1/4,1/4) sites
    assert _sites(find_interstices(nacl.copy(), tol=1e-3)) == {'8c': 8, '24d': 24}
    assert _sites(find_interstices(bulk('Si', cubic=True), tol=1e-3)) == {'8b': 8, '16c': 16, '16d': 16}
    assert _sites(find_interstices(tio2.copy(), tol=1e-3)) == {'2b': 2, '4c': 4, '4d': 4}
    # another file
    f = tmp_path / "small.dat"
    f.write_text(WYCKOFF)
    assert _sites(find_interstices(cu_bulk.copy(), tol=1e-3, wyckoff_file=str(f))) == {'4b': 4, '8c': 8}
    monkeypatch.setenv("CONFIG_ENUM_WYCKOFF", str(f))
    assert _sites(find_interstices(cu_bulk.copy(), tol=1e-3)) == {'4b': 4, '8c': 8}
    with pytest.raises(FileNotFoundError):
        find_interstices(cu_bulk.copy(), wyckoff_file=str(tmp_path / "missing.dat"))
