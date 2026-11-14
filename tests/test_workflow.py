import json
import os

import ase.io
import numpy as np
import pytest
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

from config_enum.calculators import CalcSettings, get_engine
from config_enum.workflow import EnumWorkflow, WorkflowConfig


def run(tmp_path, inpf="Cu.cif", settings=None, prog="emt", **kw):
    cfg = WorkflowConfig(base_dir=str(tmp_path), **kw)
    return EnumWorkflow(get_engine(prog, settings or CalcSettings()), cfg).run(inpf)


def test_single_step(tmp_path, cu_file):
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2])
    assert len(res.steps) == 1
    configs = res.configs
    assert len(configs) == 5 and sum(c.degeneracy for c in configs) == 120
    assert sum(c.prob for c in configs) == pytest.approx(1.)
    for i, c in enumerate(configs):
        assert c.file == os.path.join('Step1', 'input_0', 'config_%d.res' % i)
        assert ase.io.read(tmp_path / c.file).get_chemical_formula() == 'Au2Cu14'
        assert c.E_reduced == pytest.approx(c.energy - 8.6173303e-5 * 298.15 * np.log(c.degeneracy))
        # the .res file of the result has the energy in its title line
        title = (tmp_path / c.result_file).read_text().splitlines()[0].split()
        assert title[1] == 'Au2Cu14'.replace('Au2Cu14', 'AuCu7') and float(title[4]) == pytest.approx(c.energy, abs=1e-5)
    text = (tmp_path / 'Step1' / 'degeneracies').read_text()
    assert 'Config-ID   Degens      E_EMT' in text and 'Picked configs for Step 1' in text
    picked = res.picked
    assert len(picked) == 1 and picked[0].E_reduced == min(c.E_reduced for c in configs)
    assert os.listdir(tmp_path / 'Step1' / 'Picked') == ['Structures-are-not-optimised', 'input_0-' + os.path.basename(picked[0].file)] \
        or sorted(os.listdir(tmp_path / 'Step1' / 'Picked')) == sorted(['Structures-are-not-optimised', 'input_0-' + os.path.basename(picked[0].file)])
    data = json.loads((tmp_path / 'results.json').read_text())
    assert len(data['steps'][0]['configs']) == 5 and data['steps'][0]['picked'] == [picked[0].file]
    assert not (tmp_path / 'Step1' / 'Structures-are-optimised').exists()  # single points


def test_existing_output(tmp_path, cu_file):
    kw = dict(rep=(2, 1, 1), aType=['Cu'], dType='Va', noDopants=[1])
    res = run(tmp_path, cu_file, **kw)
    assert res.configs[0].atoms.get_chemical_formula() == 'Cu7'
    with pytest.raises(RuntimeError, match="use -ow"):
        run(tmp_path, cu_file, **kw)
    (tmp_path / 'Step1' / 'marker').write_text('x')
    run(tmp_path, cu_file, restart=True, **kw)
    assert (tmp_path / 'Step1' / 'marker').exists()
    run(tmp_path, cu_file, overwrite=True, **kw)
    assert not (tmp_path / 'Step1' / 'marker').exists()


def test_several_numbers_of_dopants(tmp_path, cu_file):
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[0, 1, 2])
    by_nD = {n: [c for c in res.configs if c.nD == n] for n in (0, 1, 2)}
    assert [len(by_nD[n]) for n in (0, 1, 2)] == [1, 1, 5]
    assert by_nD[0][0].atoms.get_chemical_formula() == 'Cu16' and by_nD[1][0].degeneracy == 16
    for n in (0, 1, 2):  # the statistics is done for each composition separately
        assert sum(c.prob for c in by_nD[n]) == pytest.approx(1.)
    assert [c.name for c in res.configs] == ['config_%d' % i for i in range(7)]
    with pytest.raises(ValueError):
        run(tmp_path, cu_file, overwrite=True, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[1, 2], seq=1)


@pytest.mark.parametrize("dType", ['Au', 'Va'])
def test_sequential_doping(tmp_path, cu_file, dType):
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType=dType, noDopants=[3], seq=1, seqpick=2)
    assert [s.nDops for s in res.steps] == [[1], [1], [1]]
    assert [len(s.inputs) for s in res.steps] == [1, 1, 2]
    formula = {'Au': ['AuCu15', 'Au2Cu14', 'Au3Cu13'], 'Va': ['Cu15', 'Cu14', 'Cu13']}[dType]
    comp = SymmetryEquivalenceCheck()
    for s, f in zip(res.steps, formula):
        assert all(c.atoms.get_chemical_formula() == f for c in s.configs)
        for find in set(c.input_index for c in s.configs):
            group = [c for c in s.configs if c.input_index == find]
            # every remaining site is used exactly once
            assert sum(c.degeneracy for c in group) == 16 - (s.step - 1)
            for i in range(len(group)):
                for j in range(i):
                    assert not comp.compare(group[i].atoms, group[j].atoms)
        assert len(s.picked) == min(2, len(s.configs))
        if len(s.picked) == 2:
            assert not comp.compare(s.picked[0].atoms, s.picked[1].atoms)
    # the configurations of Step 2 are the symmetrically-unique pairs
    assert len(res.steps[1].configs) == 5
    # the sites still available are passed on correctly
    last = res.steps[-1].configs[0]
    assert len(last.atoms_list) == 13 and all(last.atoms[i].symbol == 'Cu' for i in last.atoms_list)
    assert (tmp_path / 'Step3' / 'input_1').is_dir()


def test_sequential_steps_of_two(tmp_path, cu_file):
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[3], seq=2)
    assert [s.nDops for s in res.steps] == [[2], [1]]
    assert res.steps[-1].configs[0].atoms.get_chemical_formula() == 'Au3Cu13'


def test_seqsymm_ref_merges_more(tmp_path, cu_file):
    kw = dict(rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2], seq=1)
    cur = run(tmp_path, cu_file, **kw)
    ref = run(tmp_path, cu_file, overwrite=True, seqsymm='ref', **kw)
    assert len(cur.steps[1].configs) == 5 and len(ref.steps[1].configs) == 1


@pytest.mark.parametrize("seqtype", ['E', 'S', 'D', 'R'])
def test_sorting(tmp_path, cu_file, seqtype):
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2], seqtype=seqtype, seqpick=2)
    configs, picked = res.configs, res.picked
    assert len(picked) == 2
    if seqtype == 'E':
        assert [c.E_reduced for c in picked] == sorted(c.E_reduced for c in configs)[:2]
    elif seqtype == 'D':
        assert picked[0].degeneracy == max(c.degeneracy for c in configs)
    elif seqtype == 'S':
        from config_enum.symmetry import spacegroup_number
        assert spacegroup_number(picked[0].spacegroup) == max(spacegroup_number(c.spacegroup) for c in configs)


def test_output_options(tmp_path, cu_file):
    res = run(tmp_path, cu_file, rep=(2, 1, 1), aType=['Cu'], dType='Au', noDopants=[1], otype='vasp', sepfol=True)
    assert res.configs[0].file == os.path.join('Step1', 'input_0', '0', 'config_0.vasp')
    assert 'Direct' in (tmp_path / res.configs[0].file).read_text()
    res = run(tmp_path, cu_file, overwrite=True, rep=(2, 1, 1), aType=['Cu'], dType='Au', noDopants=[1], otype='vasp', sepfol=True, outf='POSCAR')
    assert res.configs[0].file == os.path.join('Step1', 'input_0', '0', 'POSCAR')
    res = run(tmp_path, cu_file, overwrite=True, rep=(2, 1, 1), aType=['Cu'], dType='Au', noDopants=[1], otype='xyz', outf='struct_x')
    assert res.configs[0].file == os.path.join('Step1', 'input_0', 'struct_0.xyz')


def test_copy_inputs_and_casm(tmp_path, cu_file):
    (tmp_path / 'INPUTS').mkdir()
    (tmp_path / 'INPUTS' / 'INCAR').write_text('ENCUT = 400\n')
    res = run(tmp_path, cu_file, rep=(2, 2, 1), aType=['Cu'], dType='Va', noDopants=[2], casm=True, copyinp=True)
    d = tmp_path / 'Step1' / 'input_0' / '1'
    assert (d / 'POS').read_text().splitlines()[0].strip() == 'config_1'
    assert (d / 'calctype.default' / 'INCAR').exists() and (d / 'calctype.default' / 'POSCAR').exists()
    prim = json.loads((tmp_path / 'Step1' / 'prim.json').read_text())
    assert len(prim['basis']) == 16 and prim['basis'][0]['occupant_dof'] == ['Cu', 'Va'] and prim['coordinate_mode'] == 'Fractional'
    cl = json.loads((tmp_path / 'Step1' / 'config_list.json').read_text())['supercells']['SCEL1_1_1_1_0_0_0']
    assert len(cl) == len(res.configs) == 5
    assert all(sum(c['dof']['occupation']) == 2 for c in cl.values())


def test_selection_options(tmp_path, nacl):
    ase.io.write(tmp_path / 'NaCl.cif', nacl)
    # Wyckoff sites, primitive cell with the labels of the conventional cell, supercell
    res = run(tmp_path, 'NaCl.cif', prog='lj', ifPrim=True, wycklist=['4a'], dType='K', noDopants=[1])
    assert len(res.configs) == 1 and res.configs[0].degeneracy == 1 and res.configs[0].atoms.get_chemical_formula() == 'ClK'
    res = run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, ifPrim=True, rep=(2, 2, 2), aType=['Na'], dType='K', noDopants=[1])
    assert len(res.configs) == 1 and res.configs[0].degeneracy == 8 and res.configs[0].atoms.get_chemical_formula() == 'Cl8KNa7'
    # atom list and excluded atoms that are doped in all configs
    atoms = ase.io.read(tmp_path / 'NaCl.cif')
    na = [a.index + 1 for a in atoms if a.symbol == 'Na']
    res = run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, list=",".join(map(str, na)), excludelist=[str(na[0])], dopeexcluded=True,
              dexcl_type='Li', dType='K', noDopants=[1])
    assert all(c.atoms.get_chemical_formula() == 'Cl4KLiNa2' for c in res.configs) and sum(c.degeneracy for c in res.configs) == 3
    # excluded atoms are left alone without -dexcl
    res = run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, aType=['Na'], excludelist=[str(na[0])], dType='K', noDopants=[1])
    assert all(c.atoms.get_chemical_formula() == 'Cl4KNa3' for c in res.configs)
    with pytest.raises(RuntimeError, match="No atoms were selected"):
        run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, aType=['Fe'], dType='K', noDopants=[1])
    # no symmetry
    res = run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, nosymm=True, aType=['Na'], dType='K', noDopants=[1])
    assert len(res.configs) == 4
    # final check with ASE
    res = run(tmp_path, 'NaCl.cif', prog='lj', overwrite=True, nosymm=True, ASE_symm=True, aType=['Na'], dType='K', noDopants=[1])
    assert len(res.configs) == 1 and res.configs[0].degeneracy == 4


def test_reference_structure(tmp_path, nacl):
    """A (partially) doped structure with the pristine structure as the reference for the sites and symmetry."""
    ase.io.write(tmp_path / 'ref.cif', nacl)
    doped = ase.io.read(tmp_path / 'ref.cif')
    doped[[a.index for a in doped if a.symbol == 'Na'][0]].symbol = 'K'
    ase.io.write(tmp_path / 'doped.cif', doped)
    res = run(tmp_path, 'doped.cif', prog='lj', aType=['Na'], wycklist=['4a'], initCrysFile='ref.cif', dType='Li', noDopants=[1])
    assert len(res.configs) == 1 and res.configs[0].degeneracy == 3


def test_dry_run(tmp_path, cu_file, capsys):
    res = run(tmp_path, cu_file, dType='Au', noDopants=[1], dry=True)
    assert res.steps == [] and not (tmp_path / 'Step1').exists()
    out = capsys.readouterr().out
    assert "4a: ['Cu0', 'Cu1', 'Cu2', 'Cu3']" in out and 'Dry run' in out


def test_optimisation_and_parallel_tasks(tmp_path, cu_file):
    s = CalcSettings(optimise=True, fmax=0.05, steps=50)
    res = run(tmp_path, cu_file, settings=s, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2], ntasks=2)
    serial = run(tmp_path, cu_file, settings=s, overwrite=True, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2])
    assert (tmp_path / 'Step1' / 'Structures-are-optimised').exists()
    assert [c.energy for c in res.configs] == pytest.approx([c.energy for c in serial.configs], abs=1e-6)
    single = run(tmp_path, cu_file, overwrite=True, rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2])
    assert all(o.energy < sp.energy for o, sp in zip(res.configs, single.configs))
    c = serial.configs[0]
    assert not np.allclose(c.final_atoms.positions, c.atoms.positions)  # relaxed result, unrelaxed structure passed on


def test_failed_calculations(tmp_path, cu_file):
    from ase.calculators.emt import EMT

    from config_enum.calculators import AseEngine

    class Flaky(AseEngine):
        n = 0

        def run(self, atoms, name="config", directory="."):
            Flaky.n += 1
            if Flaky.n == 2:
                raise RuntimeError("SCF not converged")
            return super().run(atoms, name, directory)
    kw = dict(rep=(2, 2, 1), aType=['Cu'], dType='Au', noDopants=[2], base_dir=str(tmp_path))
    res = EnumWorkflow(Flaky(EMT, label='flaky'), WorkflowConfig(**kw)).run(cu_file)
    assert [c.failed for c in res.configs] == [False, True, False, False, False]
    assert res.configs[1].energy == 0.0 and res.configs[1].spacegroup == 'P1 (1)' and not res.picked[0].failed
    Flaky.n = 0
    with pytest.raises(RuntimeError, match="stopping as requested"):
        EnumWorkflow(Flaky(EMT, label='flaky'), WorkflowConfig(overwrite=True, noskip=True, **kw)).run(cu_file)


def test_interstitials(tmp_path, cu_file):
    """-ai -at X: H on the interstices of fcc Cu (octahedral, tetrahedral and 24d sites), then two H by sequential doping."""
    res = run(tmp_path, cu_file, prog='lj', aType=['X'], addint=True, dType='H', noDopants=[1])
    assert sorted(c.degeneracy for c in res.configs) == [4, 8, 24]
    for c in res.configs:
        assert sorted(c.atoms.get_chemical_symbols()) == ['Cu'] * 4 + ['H']  # the unused interstices are removed
    res = run(tmp_path, cu_file, prog='lj', overwrite=True, aType=['X'], wycklist=['4b'], addint=True, dType='H', noDopants=[1])
    assert [c.degeneracy for c in res.configs] == [4]
    # sequential doping: the interstices that are left are available in the next step
    res = run(tmp_path, cu_file, prog='lj', overwrite=True, aType=['X'], wycklist=['4b'], addint=True, dType='H', noDopants=[2], seq=1)
    assert len(res.steps) == 2 and sum(c.degeneracy for c in res.steps[1].configs) == 3
    assert all(sorted(c.atoms.get_chemical_symbols()) == ['Cu'] * 4 + ['H'] * 2 for c in res.steps[1].configs)
    with pytest.raises(FileNotFoundError):
        run(tmp_path, cu_file, prog='lj', overwrite=True, aType=['X'], addint=True, wyckoff_file=str(tmp_path / "none.dat"), dType='H', noDopants=[1])


def test_several_input_files(tmp_path, cu_bulk):
    ase.io.write(tmp_path / "a.cif", cu_bulk * (2, 1, 1))
    ase.io.write(tmp_path / "b.cif", cu_bulk * (1, 2, 1))  # the same structure, rotated
    res = run(tmp_path, ["a.cif", "b.cif"], aType=['Cu'], dType='Au', noDopants=[2], seq=1, seqpick=3)
    assert [c.input_index for c in res.steps[0].configs] == [0, 1]
    assert len(res.steps[0].picked) == 2 and (tmp_path / 'Step2' / 'input_1').is_dir()
    assert {c.root for c in res.steps[1].configs} == {0, 1}
    text = (tmp_path / 'Step2' / 'degeneracies').read_text()
    assert 'Step1/input_0/config_0.res' in text and 'Step1/input_1/config_0.res' in text


def test_partial_occupancies(tmp_path, capsys):
    from ase.spacegroup import crystal

    from config_enum.io import partial_occupancies, read_struct
    atoms = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225, cellpar=[5.64] * 3 + [90] * 3, occupancies=[0.5, 1])
    ase.io.write(tmp_path / "NaCl.cif", atoms)
    assert partial_occupancies(read_struct(str(tmp_path / "NaCl.cif"))) == [('Na', 0.5, 4)]
    res = run(tmp_path, 'NaCl.cif', prog='lj', aType=['Na'], dType='Va', noDopants=[2])
    assert "Na: 4 sites with occupancy 0.5 (2.0 atoms; 2.0 to remove)" in capsys.readouterr().out
    assert len(res.configs) == 1 and res.configs[0].degeneracy == 6 and res.configs[0].atoms.get_chemical_formula() == 'Cl4Na2'
