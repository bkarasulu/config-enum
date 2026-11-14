"""Run the tutorials (the MACE one only when mace-torch and its model are available and CONFIG_ENUM_SLOW_TESTS=1)."""
import json
import os
import shutil
import subprocess
import sys

import pytest
from conftest import ROOT

TUT = os.path.join(ROOT, "Tutorials")


def _copy(name, tmp_path):
    dst = tmp_path / name
    shutil.copytree(os.path.join(TUT, name), dst, ignore=shutil.ignore_patterns("Step*", "results.json", "*.log"))
    return dst


def _config_enum(args, cwd, timeout=1800):
    cmd = [sys.executable, os.path.join(ROOT, "config_enum.py")] + args
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r.stdout


def test_tutorial_01_cu_au_emt(tmp_path):
    d = _copy("01_Cu-Au_EMT", tmp_path)
    _config_enum(["--input", "input.toml"], d)
    steps = json.loads((d / "results.json").read_text())["steps"]
    assert [len(s["configs"]) for s in steps] == [1, 5, 13, 24]
    step2 = steps[1]["configs"]
    assert [c["degeneracy"] for c in step2] == [12, 3, 12, 3, 1]
    assert [c["energy"] for c in step2] == pytest.approx([-0.10295, -0.19826, -0.16373, -0.13619, -0.16204], abs=2e-4)
    assert steps[1]["picked"] == ["Step2/input_0/config_2.res", "Step2/input_0/config_1.res"]
    assert all(len(s["picked"]) == 2 for s in steps[1:])
    # no sequential doping: the same pairs, with the degeneracies of the pairs
    _config_enum(["--input", "input.toml", "-ow", "-seq", "0", "-nd", "2"], d)
    configs = json.loads((d / "results.json").read_text())["steps"][0]["configs"]
    assert [c["degeneracy"] for c in configs] == [192, 48, 192, 48, 16]
    assert [c["energy"] for c in configs] == pytest.approx([c["energy"] for c in step2], abs=2e-4)


def test_tutorial_03_python_api(tmp_path):
    d = _copy("03_python_api", tmp_path)
    r = subprocess.run([sys.executable, "python_api.py"], cwd=d, capture_output=True, text=True, timeout=600,
                       env=dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src")))
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    assert "5 symmetrically-unique pairs, degeneracies [192, 48, 192, 48, 16]" in r.stdout
    assert "Step 3: 13 configurations" in r.stdout


def test_tutorial_02_sites_and_counts(tmp_path):
    """The enumeration of Tutorial 2 (with the LJ potential in place of MACE)."""
    d = _copy("02_LiCoO2_vacancies_MACE", tmp_path)
    out = _config_enum(["-i", "LiCoO2.cif", "-r", "2", "2", "1", "-dry", "-dt", "Va", "-nd", "1"], d)
    assert "12b: ['Li12'" in out and "Space group= R-3m (166)" in out
    _config_enum(["-i", "LiCoO2.cif", "-r", "2", "2", "1", "-at", "Li", "-wl", "12b", "-dt", "Va", "-nd", "1:3", "-pot", "LJ"], d)
    configs = json.loads((d / "results.json").read_text())["steps"][0]["configs"]
    assert [c["degeneracy"] for c in configs] == [12, 36, 18, 12, 24, 72, 72, 36, 12, 4]


@pytest.mark.skipif(not os.environ.get("CONFIG_ENUM_SLOW_TESTS"), reason="set CONFIG_ENUM_SLOW_TESTS=1 to run the MACE tutorial")
def test_tutorial_02_mace(tmp_path):
    pytest.importorskip("mace")
    d = _copy("02_LiCoO2_vacancies_MACE", tmp_path)
    _config_enum(["--input", "input.toml"], d, timeout=3600)
    configs = json.loads((d / "results.json").read_text())["steps"][0]["configs"]
    assert [c["degeneracy"] for c in configs] == [12, 36, 18, 12, 24, 72, 72, 36, 12, 4]
    assert all(c["energy"] < 0 and not c["failed"] for c in configs)
    for n in (1, 2, 3):
        assert sum(c["prob"] for c in configs if c["nD"] == n) == pytest.approx(1.)


def _configs(d, step=0):
    return json.loads((d / "results.json").read_text())["steps"][step]["configs"]


def test_tutorial_05_li2s_defects(tmp_path):
    """The three runs of Tutorial 5 (with the LJ potential in place of MACE)."""
    d = _copy("05_Li2S_defects", tmp_path)
    lj = ["-pot", "LJ", "-ow"]
    _config_enum(["-i", "Li2S.cif", "-r", "2", "2", "2", "-at", "Li", "-dt", "Va", "-nd", "2"] + lj, d)
    assert [c["degeneracy"] for c in _configs(d)] == [384, 192, 128, 96, 384, 128, 384, 96, 192, 32]
    _config_enum(["-i", "Li2S.cif", "-r", "2", "2", "2", "-at", "S", "-dt", "Cl", "-nd", "1", "-mod", "remove Li nearest 1 of Cl"] + lj, d)
    assert [(c["degeneracy"], c["formula"]) for c in _configs(d)] == [(32, "ClLi63S31")]
    out = _config_enum(["-i", "Li2S.cif", "-ai", "-at", "X", "-dt", "Li", "-nd", "1"] + lj, d)
    assert [(c["degeneracy"], c["formula"]) for c in _configs(d)] == [(4, "Li9S4"), (24, "Li9S4")]
    _config_enum(["-i", "Li2S.cif", "-ai", "-r", "2", "2", "2", "-at", "X", "-wl", "32b", "-dt", "Li", "-nd", "1",
                  "-mod", "remove Li nearest 1 of dopant"] + lj, d)
    assert [(c["degeneracy"], c["formula"]) for c in _configs(d)] == [(32, "Li64S32")]
    assert "New space group after adding all interstices:  Fm-3m (225)" in out


def test_tutorial_06_li3ps4_site_disorder(tmp_path):
    d = _copy("06_Li3PS4_site_disorder", tmp_path)
    out = _config_enum(["--input", "input.toml", "-pot", "LJ"], d)
    assert "Li: 4 sites with occupancy 0.28 (1.1 atoms; 2.9 to remove)" in out and "70 configuration(s) found" in out
    configs = _configs(d)
    assert [c["degeneracy"] for c in configs] == [1, 8, 8, 4, 4, 4, 4, 2, 4, 2, 4, 4, 4, 8, 8, 1]
    assert all(c["formula"] == "Li12P4S16" for c in configs)
    assert configs[0]["spacegroup"] == configs[-1]["spacegroup"] == "Pnma (62)"


def test_tutorial_07_llzo(tmp_path):
    d = _copy("07_LLZO_sequential", tmp_path)
    out = _config_enum(["-i", "LLZO.cif", "-prim", "-dry", "-dt", "Va", "-nd", "1"], d)
    assert "24d: ['Li68'" in out and "96h: ['Li80'" in out and "Li: 96 sites with occupancy 0.28" in out
    # the first steps of the sequential doping (LJ in place of MACE)
    _config_enum(["--input", "input.toml", "-pot", "LJ", "-nd", "3"], d)
    steps = json.loads((d / "results.json").read_text())["steps"]
    assert [len(s["configs"]) for s in steps][0] == 2 and len(steps) == 3
    assert [c["degeneracy"] for c in steps[0]["configs"]] == [12, 48]  # the 24d and 96h sites of the primitive cell
    assert steps[2]["configs"][0]["formula"] == "La12Li57O48Zr8" and len(steps[2]["picked"]) == 2
    # the example of the original script (without the CASM files): up to two vacancies on the 96h sites
    _config_enum(["-i", "LLZO.cif", "-prim", "-at", "Li", "-wl", "96h", "-dt", "Va", "-nd", "0:2", "-sf", "-ot", "vasp", "-o", "POSCAR", "-ow"], d)
    assert [sum(1 for c in _configs(d) if c["nD"] == n) for n in (0, 1, 2)] == [1, 1, 33]
    assert (d / "Step1" / "input_0" / "34" / "POSCAR").exists()
