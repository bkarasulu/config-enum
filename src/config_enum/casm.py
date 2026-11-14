"""Input files for CASM: prim.json and config_list.json."""
import json


def prim_json(atoms, atoms_list, dType):
    """The prim.json data: every selected atom can be occupied by its own type or by the dopant/vacancy type."""
    pos = atoms.get_scaled_positions()
    basis = []
    for i in range(len(atoms)):
        at = atoms[i]
        if i in atoms_list:  # for selected atoms.
            occ = [at.symbol, dType]
        else:
            occ = [at.symbol]
        basis.append({"coordinate": [float(x) for x in pos[i]], "occupant_dof": occ})
    symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
    return {"basis": basis, "coordinate_mode": "Fractional", "description": ' '.join(symbols) + ' ',
            "lattice_vectors": [[float(x) for x in v] for v in atoms.cell], "title": ''.join(symbols)}


def config_list_json(noAtoms, confs):
    """The config_list.json data; confs: list of the atom indices occupied by the dopant/vacancy for each config."""
    configs = {}
    for cnt, conf in enumerate(confs):
        occ = [0 for x in range(noAtoms)]
        for at in conf:
            occ[at] = 1
        configs[str(cnt)] = {"calctype.default": {"ref.default": {"properties": {}}}, "dof": {"occupation": occ},
                             "selected": True, "source": [{"enumerated_by": "ConfigEnumAllOccupations", "step": cnt}]}
    return {"supercells": {"SCEL1_1_1_1_0_0_0": configs}}


def write_casm(atoms, atoms_list, dType, confs, prim_file='prim.json', config_file='config_list.json'):
    with open(prim_file, 'w') as f:
        json.dump(prim_json(atoms, atoms_list, dType), f, indent=2)
    with open(config_file, 'w') as f:
        json.dump(config_list_json(len(atoms), confs), f, indent=2)
