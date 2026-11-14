"""Simple classical potentials of ASE (for quick tests): Lennard-Jones, Morse and EMT."""

ASE_POTS = ("LJ", "Morse", "EMT")


def classical_calculator(name="LJ"):
    """ASE calculator of a classical potential.

    LJ and Morse potentials (with the default ASE parameters) do not differentiate between the diff atom types, so
    give the same energies for substitutional configurations; EMT covers Al, Cu, Ag, Au, Ni, Pd, Pt (and H, C, N, O).
    VASP/CASTEP or a machine-learning potential should be used instead for real work.
    """
    key = name.strip().lower()
    if key in ("lj", "lennardjones", "lennard-jones"):
        from ase.calculators.lj import LennardJones as LJ
        return LJ()
    if key in ("morse", "mp"):
        from ase.calculators.morse import MorsePotential as MP
        return MP()
    if key == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    raise ValueError("Unknown ASE potential '%s'; choose from %s" % (name, ", ".join(ASE_POTS)))
