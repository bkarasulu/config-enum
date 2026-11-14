import numpy as np
import pytest

from config_enum.geometry import compCells, dist, lenAngleCell, lenCell, minDist, volume
from config_enum.utils import cartesian, combinations2, hashid, ncr, parse_ranges, sort_dict_by_val, sorted_nicely


def test_ncr_and_combinations():
    assert ncr(5, 2) == 10 and ncr(96, 2) == 4560 and ncr(4, 0) == 1 and ncr(4, 4) == 1 and ncr(2, 3) == 0
    assert list(combinations2('ABCD', 2)) == [('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'), ('B', 'D'), ('C', 'D')]
    assert len(list(combinations2(range(6), 3))) == ncr(6, 3)


def test_cartesian():
    out = cartesian(([1, 2, 3], [4, 5], [6, 7]))
    assert out.shape == (12, 3) and list(out[0]) == [1, 4, 6] and list(out[-1]) == [3, 5, 7]


def test_sorted_nicely_and_hashid():
    assert sorted_nicely(['24d', '8a', '96h', '16b']) == ['8a', '16b', '24d', '96h']
    assert sorted_nicely(['Li10', 'Li2', 'Li1']) == ['Li1', 'Li2', 'Li10']
    assert len(hashid()) == 6 and len(hashid(10)) == 10


def test_sort_dict_by_val():
    assert sort_dict_by_val({'a': 3., 'b': 1., 'c': 2.}) == ['b', 'c', 'a']
    assert sort_dict_by_val({'a': [0, 3.], 'b': [0, 1.]}, index=1) == ['b', 'a']
    assert sort_dict_by_val({'a': 1}) == ['a'] and sort_dict_by_val({}) == []


def test_parse_ranges():
    assert parse_ranges(['0:4', '24']) == [0, 1, 2, 3, 4, 24]
    assert parse_ranges('1:3,7') == [1, 2, 3, 7]
    with pytest.raises(ValueError):
        parse_ranges(['a'])


def test_geometry():
    cell = np.diag([2., 3., 4.])
    assert volume(cell) == pytest.approx(24.)
    assert lenCell(cell) == pytest.approx([2., 3., 4.])
    assert lenAngleCell(cell) == pytest.approx([2., 3., 4., 90., 90., 90.])
    assert compCells(2 * cell, cell) == (2, 2, 2)
    assert dist([0, 0, 0], [3, 4, 0]) == pytest.approx(5.)
    assert minDist([0.1, 0, 0], [1.9, 0, 0], cell) == pytest.approx(0.2)
