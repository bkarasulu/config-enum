import numpy as np
import pytest

from config_enum.thermo import boltz_dist, config_thermo, kb


def test_equal_energies():
    """Equal energies: the probabilities are proportional to the degeneracies and the entropy is the maximum one."""
    res = config_thermo([-1.0, -1.0, -1.0], T=300., omega=[1, 2, 5])
    assert res["P"] == pytest.approx([1 / 8, 2 / 8, 5 / 8])
    assert res["S"] == pytest.approx(kb * np.log(8)) and res["S_max"] == pytest.approx(kb * np.log(8))
    assert res["E_avg"] == pytest.approx(-1.0)
    assert res["F"] == pytest.approx(-1.0 - kb * 300. * np.log(8))


def test_against_complete_space():
    """The reduced-space statistics equals the plain Boltzmann statistics over all (degenerate) configurations."""
    E = [-10.00, -9.98, -9.95, -9.90]
    om = [2, 6, 1, 12]
    T = 500.
    res = config_thermo(E, T=T, omega=om)
    allE = np.repeat(E, om)
    w = np.exp(-(allE - allE.min()) / kb / T)
    Z = w.sum()
    p = w / Z
    assert sum(res["P"]) == pytest.approx(1.)
    assert res["P"] == pytest.approx([p[allE == e].sum() for e in E])
    assert res["E_avg"] == pytest.approx((p * allE).sum())
    assert res["F"] == pytest.approx(allE.min() - kb * T * np.log(Z))
    assert res["S"] == pytest.approx(-kb * (p * np.log(p)).sum())
    assert res["S"] <= res["S_max"]
    assert res["E_reduced"] == pytest.approx([e - kb * T * np.log(o) for e, o in zip(E, om)])


def test_boltz_dist():
    P, Em_bar = boltz_dist([0.0, 0.1], T=298.15)  # no degeneracies given: all 1
    assert Em_bar == pytest.approx([0.0, 0.1])
    assert P[1] / P[0] == pytest.approx(np.exp(-0.1 / kb / 298.15))
    P, Em_bar = boltz_dist([-5.0], omega=[4])
    assert P == pytest.approx([1.0])
    # a large energy difference must not fail
    P, _ = boltz_dist([0.0, 1000.0], T=10.)
    assert P == pytest.approx([1.0, 0.0])
    with pytest.raises(ValueError):
        boltz_dist([0., 1.], omega=[1])
