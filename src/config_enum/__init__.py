"""config_enum: automated enumeration of defect (vacancy/dopant/interstitial) configurations.

Main modules:
    config_enum.symmetry     symmetry operations, equivalent (Wyckoff) sites, interstices
    config_enum.enumeration  atom selection, symmetry-unique configurations, doping and final modifications
    config_enum.thermo       configurational thermodynamics (Boltzmann statistics with degeneracies)
    config_enum.calculators  engines for VASP, CASTEP and ASE calculators (classical and machine-learning potentials)
    config_enum.workflow     the (sequential) enumeration workflow
    config_enum.cli          the command line (``config_enum`` command)

Copyright (C) Bora Karasulu and contributors.
config_enum is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version. It is distributed WITHOUT ANY WARRANTY; see the LICENSE file for details.
"""
__version__ = "1.0.0.dev0"

from .calculators import AseEngine, CalcSettings, get_engine  # noqa: F401
from .enumeration import apply_dopants, find_unique_configs, select_atoms  # noqa: F401
from .symmetry import find_equisites, find_interstices, get_equivalent_sites, get_symmetry  # noqa: F401
from .thermo import boltz_dist, config_thermo  # noqa: F401
