import json
import os
import subprocess
import sys

import pytest

from config_enum.cli import (apply_input_file, build_parser, config_from_args, get_prog, load_input_file, main, make_settings,
                             parse_hubbard_u, parse_no_dopants)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_parse_helpers():
    assert parse_hubbard_u(['Fe', 'd', '2.7,', 'Sm', 'f', '6.1']) == {'Fe': {'d': 2.7}, 'Sm': {'f': 6.1}}
    assert parse_hubbard_u(None) == {}
    with pytest.raises(ValueError):
        parse_hubbard_u(['Fe', '2.7'])
    assert parse_no_dopants(['0:2']) == [0, 1, 2] and parse_no_dopants(['4', '1:2']) == [1, 2, 4]
    for bad in (['x'], [], ['-1']):
        with pytest.raises(ValueError):
            parse_no_dopants(bad)


def test_original_options_are_kept():
    """All options of the original script (config_enum-v2.py) are accepted."""
    args = build_parser().parse_args(
        "-i LLZO.cif -prim -at Li -wl 96h -dt Va -nd 0:2 --casm -sf -ot vasp -o POS -it cif -od OUT -ow -rs -iseed seed -ref ref.cif "
        "-tol 1e-2 -nosymm -r 2 2 1 -nl 1,2,5:8 -excl 1:4 24 -dexcl -dexcl_type Na -mod x -ai -ase -castep -vasp -vexe vasp_gam "
        "-cexe castep -copy -id INP -mp --potcarxc potpaw_PBE -xc pbesol -hubU Fe d 2.7 -mgm Fe 5.0 -np 16 -nt 2 -nn 1 -mpi -t 500 "
        "-seq 1 -spick 2 -stype S -v -dry -noskip".split())
    assert args.inpf == ['LLZO.cif'] and args.ifPrim and args.aType == ['Li'] and args.wycklist == ['96h'] and args.dType == 'Va'
    assert args.noDopants == ['0:2'] and args.rep == [2, 2, 1] and args.excludelist == ['1:4', '24'] and args.seqtype == 'S'
    assert get_prog(args) == 'vasp'
    s = make_settings(args)
    assert s.hubU == {'Fe': {'d': 2.7}} and s.magmoms == ['Fe', '5.0'] and s.idir == 'INP' and s.xc == 'pbesol' and s.mpirun
    cfg = config_from_args(args)
    assert cfg.noDopants == [0, 1, 2] and cfg.rep == (2, 2, 1) and cfg.temp == 500. and cfg.seqpick == 2 and cfg.casm and cfg.dry


def test_engine_selection():
    p = build_parser()
    base = "-i a.cif -dt Va -nd 1 -at Li".split()
    assert get_prog(p.parse_args(base)) == 'ase'
    assert get_prog(p.parse_args(base + ['-castep'])) == 'castep'
    assert get_prog(p.parse_args(base + ['--mlip', 'mace'])) == 'mlip'
    args = p.parse_args(base + ['-mlip', 'chgnet', '--mlip_kwargs', "{'a': 1}", '--relax_cell', '--device', 'cpu'])
    s = make_settings(args)
    assert s.mlip == 'chgnet' and s.mlip_kwargs == {'a': 1} and s.optimise and s.relax_cell and s.device == 'cpu'
    with pytest.raises(ValueError):
        make_settings(p.parse_args(base + ['--mlip_kwargs', "[1]"]))
    with pytest.raises(ValueError):
        make_settings(p.parse_args(base + ['-mgm', 'Fe']))


def test_main_and_input_files(tmp_path, cu_file, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert main(f"-i {cu_file} -r 2 2 1 -at Cu -dt Au -nd 2 -pot EMT --save_input saved.toml".split()) == 0
    out = capsys.readouterr().out
    assert "Total of 5 symmetrically-unique configurations were generated." in out and "E_EMT" in out
    assert "config_enum -i Cu.cif" in (tmp_path / "Step1" / "degeneracies").read_text()
    first = json.loads((tmp_path / "results.json").read_text())
    # existing output: stop with a message
    assert main(f"-i {cu_file} -r 2 2 1 -at Cu -dt Au -nd 2".split()) == 1
    assert "use -ow" in capsys.readouterr().out
    # the saved input file reproduces the run; command-line options override the file
    saved = load_input_file("saved.toml")
    assert saved == {'aType': ['Cu'], 'ase_pot': 'EMT', 'dType': 'Au', 'inpf': ['Cu.cif'], 'noDopants': ['2'], 'rep': [2, 2, 1]}
    assert main("--input saved.toml -ow -q".split()) == 0
    second = json.loads((tmp_path / "results.json").read_text())
    assert [c['energy'] for c in second['steps'][0]['configs']] == pytest.approx([c['energy'] for c in first['steps'][0]['configs']])
    assert main("--input saved.toml -ow -q -nd 1".split()) == 0
    assert len(json.loads((tmp_path / "results.json").read_text())['steps'][0]['configs']) == 1


def test_input_file_formats(tmp_path):
    (tmp_path / "in.toml").write_text('inpf = "a.cif"\n[selection]\nat = ["Li"]\ndType = "Va"\nnd = "0:2"\n[calculation]\n'
                                      'mlip = "mace"\nmlip_kwargs = {default_dtype = "float32"}\nopt = true\nseq = 1\nr = [2, 2, 1]\nnl = [1, 2, 5]\n')
    p = build_parser()
    apply_input_file(p, str(tmp_path / "in.toml"))
    args = p.parse_args([])
    assert args.inpf == ['a.cif'] and args.aType == ['Li'] and args.noDopants == ['0:2'] and args.optimise and args.seq == 1
    assert args.rep == [2, 2, 1] and args.list == '1,2,5' and make_settings(args).mlip_kwargs == {'default_dtype': 'float32'}
    (tmp_path / "in.json").write_text(json.dumps({"inpf": ["a.cif"], "dType": "Na", "noDopants": [1, 2], "at": "Li"}))
    p = build_parser()
    apply_input_file(p, str(tmp_path / "in.json"))
    assert config_from_args(p.parse_args([])).noDopants == [1, 2]
    (tmp_path / "bad.toml").write_text('unknown_option = 1\n')
    with pytest.raises(ValueError, match="Unknown option"):
        apply_input_file(build_parser(), str(tmp_path / "bad.toml"))
    (tmp_path / "bad2.toml").write_text('opt = "yes"\n')
    with pytest.raises(ValueError, match="true or false"):
        apply_input_file(build_parser(), str(tmp_path / "bad2.toml"))


def test_errors(tmp_path, cu_file, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        main(f"-i {cu_file} -dt Au -nd 1".split())  # no atom selection
    with pytest.raises(SystemExit):
        main(f"-i {cu_file} -at Cu -dt Au -nd 1 2 -seq 1".split())
    with pytest.raises(SystemExit):
        main(f"-i {cu_file} -at Cu -dt Au -nd 1 --mlip not-a-potential".split())
    capsys.readouterr()
    assert main("-i missing.cif -at Cu -dt Au -nd 1".split()) == 1
    assert "Error in reading the input file" in capsys.readouterr().out


def test_dry_run_and_script(tmp_path, cu_file):
    """The config_enum.py script in the repository takes the same options."""
    out = subprocess.run([sys.executable, os.path.join(ROOT, "config_enum.py"), "-i", cu_file, "-dry", "-dt", "Va", "-nd", "1"],
                         cwd=tmp_path, capture_output=True, text=True)
    assert out.returncode == 0 and "4a: ['Cu0', 'Cu1', 'Cu2', 'Cu3']" in out.stdout and "Space group= Fm-3m (225)" in out.stdout
    out = subprocess.run([sys.executable, os.path.join(ROOT, "config_enum.py"), "--version"], cwd=tmp_path, capture_output=True, text=True)
    assert out.stdout.startswith("config_enum ")
