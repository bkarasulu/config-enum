"""Reading and writing structures."""
import os

import ase.io
from ase.calculators.singlepoint import SinglePointCalculator


def read_struct(inpf, itype=None):
    """Read a structure with ASE (the type is determined automatically from the extension if not given)."""
    try:
        if itype:
            return ase.io.read(inpf, format=itype)
        if os.path.basename(inpf) in ("POSCAR", "CONTCAR"):
            return ase.io.read(inpf, format="vasp")
        return ase.io.read(inpf)
    except Exception as err:
        raise RuntimeError("Error in reading the input file (%s): %s" % (inpf, err))


def write_struct(fname, atoms, otype="res"):
    """Write a structure in any format of ASE (VASP files in direct coordinates, VASP5 format)."""
    d = os.path.dirname(fname)
    if d:
        os.makedirs(d, exist_ok=True)
    if any(k.startswith('config_enum') for k in atoms.arrays):  # internal labels
        atoms = atoms.copy()
        for k in [k for k in atoms.arrays if k.startswith('config_enum')]:
            del atoms.arrays[k]
    if otype == "vasp":
        ase.io.write(fname, atoms, format=otype, direct=True, vasp5=True)
    else:
        ase.io.write(fname, atoms, format=otype)


def write_res_result(fname, atoms, energy, pressure=0., spacegroup="P1", name=None):
    """Write a SHELX .res file with the energy, pressure and space group in the TITL line (AIRSS style)."""
    import ase.io.res
    out = atoms.copy()
    out.calc = SinglePointCalculator(out, energy=energy)
    out.info = {'name': name or out.get_chemical_formula('hill', empirical=True), 'pressure': pressure,
                'energy': energy, 'spacegroup': spacegroup, 'times_found': 1}
    res = ase.io.res.Res(out, name=out.info['name'], pressure=pressure, energy=energy, spacegroup=spacegroup, times_found=1)
    res.write_file(fname, write_info=False, significant_figures=6)


def partial_occupancies(atoms):
    """Partially occupied sites of a structure read from a CIF file: [(atom type, occupancy, number of sites), ...]."""
    occ = atoms.info.get('occupancy')
    kinds = atoms.arrays.get('spacegroup_kinds')
    if not occ or kinds is None:
        return []
    counts = {}
    for k in kinds:
        for sym, o in occ.get(str(int(k)), {}).items():
            if abs(float(o) - 1.) > 1e-6:
                counts[(sym, float(o))] = counts.get((sym, float(o)), 0) + 1
    return [(sym, o, n) for (sym, o), n in sorted(counts.items())]
