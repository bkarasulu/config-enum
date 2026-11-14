import os

import ase.io
import pytest
from ase.build import bulk

from config_enum.calculators import AseEngine, CalcSettings, CastepEngine, VaspEngine, get_engine
from config_enum.calculators.castep import castep_command
from config_enum.calculators.classical import classical_calculator
from config_enum.calculators.mlip import MLIPNotAvailable, available_mlips, canonical_name, mlip_calculator
from config_enum.calculators.optimise import get_optimiser, relax
from config_enum.calculators.vasp import (copy_vasp_inputs, get_potcar, handle_magmoms, make_potcar, prepare_vasp, read_pressure,
                                          vasp_command, vasp_continue)


def emt_factory(**kwargs):
    """A 'custom' machine-learning potential for the tests: custom:test_calculators:emt_factory"""
    from ase.calculators.emt import EMT
    emt_factory.kwargs = kwargs
    return EMT()


def test_classical_calculators():
    for name in ("LJ", "Morse", "EMT"):
        atoms = bulk("Cu", cubic=True)
        atoms.calc = classical_calculator(name)
        assert isinstance(atoms.get_potential_energy(), float)
    with pytest.raises(ValueError):
        classical_calculator("ReaxFF")


def test_ase_engine_single_point_and_relax():
    atoms = bulk("Cu", cubic=True) * (2, 1, 1)
    atoms[0].symbol = 'Au'
    E0, a0, P0 = get_engine("emt").run(atoms.copy())
    E1, a1, P1 = get_engine("emt", CalcSettings(optimise=True)).run(atoms.copy())
    E2, a2, P2 = get_engine("ase", CalcSettings(ase_pot="EMT", optimise=True, relax_cell=True, optimiser="FIRE")).run(atoms.copy())
    assert E2 < E1 < E0 and abs(P2) < abs(P0)
    assert a2.get_volume() > a0.get_volume()
    assert get_engine("emt").describe() == "EMT (single point)" and "cell relaxed" in get_engine("emt", CalcSettings(optimise=True, relax_cell=True)).describe()
    # an ASE calculator (instance or class) can be used directly
    from ase.calculators.emt import EMT
    assert AseEngine(EMT()).run(atoms.copy())[0] == pytest.approx(E0)
    assert AseEngine(EMT, label="EMT").run(atoms.copy())[0] == pytest.approx(E0)
    with pytest.raises(ValueError):
        get_engine("gaussian")
    with pytest.raises(ValueError):
        get_engine("mlip")  # no potential given


def test_optimiser():
    assert get_optimiser("fire").__name__ == "FIRE"
    with pytest.raises(ValueError):
        get_optimiser("newton")
    atoms = bulk("Cu", cubic=True)
    atoms.rattle(0.05, seed=1)
    E, atoms = relax(atoms, classical_calculator("EMT"), fmax=0.01)
    assert abs(atoms.get_forces()).max() < 0.01


def test_mlip_backends():
    assert "mace" in available_mlips() and "chgnet" in available_mlips()
    assert canonical_name("MACE-MP") == "mace" and canonical_name("uma") == "fairchem"
    with pytest.raises(ValueError, match="Unknown ML potential"):
        mlip_calculator("not-a-potential")
    try:
        import mace  # noqa: F401
    except ImportError:
        with pytest.raises(MLIPNotAvailable, match="pip install mace-torch"):
            mlip_calculator("mace")


def test_custom_mlip_engine():
    s = CalcSettings(mlip="custom:test_calculators:emt_factory", mlip_model="my.model", mlip_kwargs={"scale": 2})
    engine = get_engine("mlip", s)
    atoms = bulk("Cu", cubic=True)
    E, _, _ = engine.run(atoms.copy())
    assert E == pytest.approx(get_engine("emt").run(atoms.copy())[0])
    assert emt_factory.kwargs == {"scale": 2, "model": "my.model"}
    assert engine.get_calculator() is engine.get_calculator()  # the model is loaded only once


def test_commands(monkeypatch):
    monkeypatch.delenv("VASP_COMMAND", raising=False)
    monkeypatch.delenv("CASTEP_COMMAND", raising=False)
    assert vasp_command(nprocs=32, ntasks=4, nnodes=2) == "srun -n 8 -N 1 vasp_std"
    assert vasp_command(nprocs=32, ntasks=2, mpirun=True, vaspexe="vasp_gam") == "mpirun -np 16 vasp_gam"
    assert castep_command(nprocs=16, mpirun=True) == "mpirun -np 16 castep.mpi"
    assert VaspEngine(CalcSettings(nprocs=8)).exe == "srun -n 8 -N 1 vasp_std"
    assert CastepEngine(CalcSettings(nprocs=8, mpirun=True)).cascmd == "mpirun -np 8 castep.mpi"
    monkeypatch.setenv("VASP_COMMAND", "mpirun vasp")
    monkeypatch.setenv("CASTEP_COMMAND", "castep.serial")
    assert vasp_command() == "mpirun vasp" and castep_command() == "castep.serial"
    assert vasp_command("my_vasp") == "my_vasp"


def test_magmoms():
    atoms = bulk("Cu", cubic=True)
    atoms[0].symbol = 'Fe'
    handle_magmoms(atoms, ['Fe', '5', 'O', '0.6'])
    assert atoms.get_initial_magnetic_moments().tolist() == [5., 0., 0., 0.]
    assert handle_magmoms(atoms, None) is atoms


def test_potcar(tmp_path, fake_potcars, monkeypatch):
    atoms = bulk("Cu", cubic=True)
    atoms[0].symbol = 'Au'
    ase.io.write(tmp_path / "POSCAR", atoms, format="vasp")
    make_potcar(wDir=str(tmp_path))
    text = (tmp_path / "POTCAR").read_text()
    assert text.count("TITEL") == 2 and text.index("PAW_PBE Au") < text.index("PAW_PBE Cu")
    assert get_potcar(["Cu"], "potpaw_PBE").count("TITEL") == 1
    monkeypatch.delenv("VASP_PP_PATH")
    with pytest.raises(RuntimeError, match="VASP_PP_PATH"):
        make_potcar(wDir=str(tmp_path))


def test_prepare_vasp(tmp_path, fake_potcars):
    """INCAR and KPOINTS of the input folder are used, and the magnetic moments and Hubbard U are added."""
    idir = tmp_path / "INPUTS"
    idir.mkdir()
    (idir / "INCAR").write_text("ENCUT = 450\nIBRION = 2\nNSW = 50\nISPIN = 2\n")
    (idir / "KPOINTS").write_text("K-Points\n 0\nMonkhorst-Pack\n 3 3 2\n 0 0 0\n")
    wdir = tmp_path / "config_0"
    wdir.mkdir()
    atoms = bulk("Cu", cubic=True)
    atoms[0].symbol = 'Fe'
    atoms[1].symbol = 'O'
    ase.io.write(wdir / "POSCAR", atoms, format="vasp")
    copy_vasp_inputs(str(idir), str(wdir))
    atoms = prepare_vasp(str(wdir), xc='pbe', mgms=['Fe', 5.0], hubU={'Fe': {'d': 2.7}})
    atoms.calc.write_input(atoms)
    incar = (wdir / "INCAR").read_text()
    assert "ENCUT = 450" in incar and "NSW = 50" in incar and "LDAU = .TRUE." in incar and "LASPH = .TRUE." in incar
    order = list(dict.fromkeys(ase.io.read(wdir / "POSCAR").get_chemical_symbols()))
    ldauu = [x for x in incar.splitlines() if "LDAUU" in x][0].split("=")[1].split()
    ldaul = [x for x in incar.splitlines() if "LDAUL" in x][0].split("=")[1].split()
    assert float(ldauu[order.index('Fe')]) == 2.7 and int(ldaul[order.index('Fe')]) == 2
    assert sum(float(x) for x in ldauu) == 2.7
    assert "MAGMOM" in incar and "5.0" in [x for x in incar.splitlines() if "MAGMOM" in x][0]
    assert "3 3 2" in (wdir / "KPOINTS").read_text()
    assert (wdir / "POTCAR").read_text().count("TITEL") == 3


def test_vasp_continue_and_pressure(tmp_path):
    for fn in ("INCAR", "POSCAR", "OUTCAR", "WAVECAR"):
        (tmp_path / fn).write_text(fn)
    (tmp_path / "CONTCAR").write_text("relaxed")
    vasp_continue(wdir=str(tmp_path))
    vasp_continue(ifDel=1, wdir=str(tmp_path))
    assert sorted(os.listdir(tmp_path / "RUN1")) == ["CONTCAR", "INCAR", "OUTCAR", "POSCAR", "WAVECAR"]
    assert (tmp_path / "RUN2").is_dir() and (tmp_path / "POSCAR").read_text() == "relaxed" and not (tmp_path / "WAVECAR").exists()
    (tmp_path / "OUTCAR").write_text("  external pressure =        1.50 kB  Pullay stress =        0.00 kB\n"
                                     "  external pressure =       -3.21 kB  Pullay stress =        1.21 kB\n")
    assert read_pressure(str(tmp_path / "OUTCAR")) == pytest.approx(-0.2)
    assert read_pressure(str(tmp_path / "missing")) == 0.
