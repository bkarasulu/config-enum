"""Configurational thermodynamics of a set of symmetry-inequivalent configurations.

Configurational statistics as given in R. Grau-Crespo et al., J. Phys.: Condens. Matter 19 (2007) 256201 and
DOI:10.13140/RG.2.1.3561.4161 (Mol2Net 2015 paper).
"""
import math

import numpy as np

from .log import debug, say

kb = 8.6173303e-5  # Boltzmann constant (eV/K).


def config_thermo(energies, T=298.15, omega=()):
    """Configurational statistics of M inequivalent configurations with energies E_m and degeneracies Omega_m.

    Higher degeneracy for an irreducible configuration means higher configurational entropy, and thus lower reduced
    energy: E_m_bar = E_m - T*S_m with S_m = kb*ln(Omega_m), and P_m = exp(-E_m_bar/kb/T)/Z.

    Returns a dict with:
        P (list): probabilities of the inequivalent configurations (degeneracies included); sum to 1.
        P_single (list): probability of ONE of the Omega_m equivalent configurations (P/Omega).
        E_reduced (list): reduced energies E_m - kb*T*ln(Omega_m) (eV).
        E_avg: average energy in configurational equilibrium (eV).
        F: configurational free energy -kb*T*ln(Z), on the scale of the input energies (eV).
        S: configurational entropy (E_avg - F)/T (eV/K).
        S_max: upper limit of the entropy, kb*ln(sum(Omega)) (eV/K).
        Z: partition function relative to the lowest reduced energy.
    """
    energies = [float(E) for E in energies]
    if len(omega) == 0:  # If the degeneracies are not given explicitly, all set to 1.
        omega = [1 for E in energies]
    omega = [float(om) for om in omega]
    if len(omega) != len(energies):
        raise ValueError("boltz_dist: %d energies but %d degeneracies" % (len(energies), len(omega)))

    Sm = [kb * T * np.log(om) for om in omega]  # degeneracy entropy (times T, eV)
    # Reduced energy; there is a typo in the Grau-Crespo, J. Phys.: Condens. Matter 19 (2007) 256201 paper, Eqn 12.
    Em_bar = [energies[i] - Sm[i] for i in range(len(energies))]
    mn = min(Em_bar)
    rel_Em_bar = [i - mn for i in Em_bar]
    Pm = []
    for E in rel_Em_bar:
        try:
            Pm.append(math.exp(-E / kb / T))
        except (OverflowError, ZeroDivisionError):
            Pm.append(0.0)
    Z_bar = sum(Pm)  # i.e. partition fnc
    Pm_bar = [P / Z_bar for P in Pm]  # reduced probability for an independent config.

    E_avg = sum(energies[i] * Pm_bar[i] for i in range(len(energies)))
    F = mn - kb * T * np.log(Z_bar)
    S = (E_avg - F) / T
    Smax = kb * np.log(sum(omega))  # should be total no of possible configs (so sum(degens))
    return {"P": Pm_bar, "P_single": [Pm_bar[i] / omega[i] for i in range(len(omega))], "E_reduced": Em_bar,
            "E_avg": E_avg, "F": float(F), "S": float(S), "S_max": float(Smax), "Z": Z_bar, "T": T,
            "T_Sm": [float(s) for s in Sm]}


def boltz_dist(energies, T=298.15, omega=(), verb=0):
    """Return the occupation probabilities of the configurations at a given temperature based on their energies,
    and their reduced energies: (Pm_bar, Em_bar)."""
    res = config_thermo(energies, T=T, omega=omega)
    say('\nConfigurational thermodynamics at T= %d K\n' % T)
    out = say if verb else debug
    out("Average energy of the system in configurational equilibrium,  E=%.5f eV" % res["E_avg"])
    out("Configurational free energy, F=%.5f eV" % res["F"])
    out("Configurational entropy, S=%.3e eV/K" % res["S"])
    out("Upper limit of config. entropy, Smax= %.3e eV/K" % res["S_max"])
    out("Degeneracy entropy (T*Sm, eV)= %s" % (['%.3f' % m for m in res["T_Sm"]]))
    out("Probabilities of a single (symmetry-equivalent) configuration: %s" % (['%.3f' % m for m in res["P_single"]]))
    out("Reduced probabilities for independent configurations: %s" % (['%.3f' % m for m in res["P"]]))
    return res["P"], res["E_reduced"]
