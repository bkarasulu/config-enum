"""Tutorial 3: config_enum from Python (EMT; replace the engine by an ML potential or VASP for real work)."""
import ase.io

from config_enum.calculators import CalcSettings, get_engine
from config_enum.enumeration import apply_dopants, enumerate_configs, select_atoms
from config_enum.symmetry import get_equivalent_sites, get_symmetry
from config_enum.thermo import config_thermo
from config_enum.workflow import EnumWorkflow, WorkflowConfig

# --- 1. the building blocks ------------------------------------------------------------------------------------
atoms = ase.io.read("Cu.cif").repeat((2, 2, 2))

data = get_symmetry(atoms, tol=1e-3)                  # symmetry operations, equivalent atoms and Wyckoff letters (spglib)
eSites, wyckoffs = get_equivalent_sites(atoms, data)
print("Space group:", data["spacegroup"], " sites:", {k: len(v) for k, v in eSites.items()})

atoms_list = select_atoms(atoms, aType=["Cu"])        # the sites that can be doped (as -at Cu)
confs, degens = enumerate_configs(atoms, atoms_list, 2, data["rotations"], data["translations"])
print("%d symmetrically-unique pairs, degeneracies %s" % (len(confs), degens))

engine = get_engine("emt", CalcSettings(optimise=True, fmax=0.02))
# engine = get_engine("mlip", CalcSettings(mlip="mace", mlip_model="small", optimise=True))
energies = []
for conf in confs:
    doped, sites_left = apply_dopants(atoms, conf, "Au", atoms_list=atoms_list)
    E, doped, P = engine.run(doped)
    energies.append(E)

thermo = config_thermo(energies, T=300., omega=degens)
for conf, E, om, p in zip(confs, energies, degens, thermo["P"]):
    print("sites %-10s E = %9.5f eV  degeneracy %3d  probability %.3f" % (conf, E, om, p))
print("Configurational entropy at 300 K: %.3e eV/K (maximum %.3e eV/K)" % (thermo["S"], thermo["S_max"]))

# --- 2. the whole workflow (as the config_enum command) ---------------------------------------------------------
cfg = WorkflowConfig(rep=(2, 2, 2), aType=["Cu"], dType="Au", noDopants=[3], seq=1, seqpick=2, overwrite=True, temp=300.)
results = EnumWorkflow(engine, cfg).run("Cu.cif")
for step in results.steps:
    best = min(step.configs, key=lambda c: c.E_reduced)
    print("Step %d: %d configurations; lowest reduced energy: %s (%s), E = %.5f eV"
          % (step.step, len(step.configs), best.file, best.spacegroup, best.energy))
print("Picked for the next step:", [c.file for c in results.picked])
