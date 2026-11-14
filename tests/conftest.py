import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(ROOT, "src"), ROOT, os.path.join(ROOT, "tests")):  # works without 'pip install -e .' as well
    if p not in sys.path:
        sys.path.insert(0, p)


@pytest.fixture
def cu_bulk():
    from ase.build import bulk
    return bulk("Cu", "fcc", a=3.61, cubic=True)


@pytest.fixture
def nacl():
    from ase.spacegroup import crystal
    return crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225, cellpar=[5.64] * 3 + [90] * 3)


@pytest.fixture
def tio2():
    from ase.spacegroup import crystal
    return crystal(['Ti', 'O'], [(0, 0, 0), (0.305, 0.305, 0)], spacegroup=136, cellpar=[4.59, 4.59, 2.96, 90, 90, 90])


@pytest.fixture
def cu_file(tmp_path, cu_bulk):
    import ase.io
    ase.io.write(tmp_path / "Cu.cif", cu_bulk)
    return "Cu.cif"


@pytest.fixture
def fake_potcars(tmp_path, monkeypatch):
    """Minimal placeholder POTCAR files, enough for ASE to write VASP inputs (no VASP run)."""
    pp = tmp_path / "pp"
    for el in ("Cu", "Au", "Fe", "O", "Li"):
        for sub in (el, el + "_pv", el + "_sv"):
            d = pp / "potpaw_PBE" / sub
            d.mkdir(parents=True)
            (d / "POTCAR").write_text(
                f"  PAW_PBE {el} 06Sep2000\n 11.0\n   VRHFIN ={el}: d s\n   TITEL  = PAW_PBE {el} 06Sep2000\n"
                "   POMASS =   63.546; ZVAL   =   11.000    mass and valenz\n   ENMAX  =  295.446; ENMIN  =  221.585 eV\n End of Dataset\n")
    monkeypatch.setenv("VASP_PP_PATH", str(pp))
    return pp
