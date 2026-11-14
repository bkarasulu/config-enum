"""The symmetry-equivalence check of the original script (config_enum-v2.py), kept as a reference for the tests.

It compares the scaled coordinates of each new configuration with the symmetry images of the irreducible
configurations found so far. `shift` reproduces the atom-index offset of the original script (scpos[at-1], shift=1);
shift=0 uses the positions of the selected atoms themselves.
"""
import numpy as np


def legacy_unique_configs(c, scpos, all_W, all_w, tol=1e-3, shift=0):
    scpos = np.where(np.isclose(scpos, 1.0), 0.0, scpos)
    all_eSites = []
    degens = []
    unique = []
    for cind, conf in enumerate(c):
        coords = np.zeros((4, len(conf)))
        for i, at in enumerate(conf):
            coords[0:3, i] = scpos[at - shift][0:3]
            coords[3, i] = 1.
        eFlag = 0
        for i, equi in enumerate(all_eSites):
            for crd in equi:
                ifFound = np.zeros((len(conf)))
                used = np.zeros((len(conf)))
                for k in range(len(conf)):
                    for j in range(len(conf)):
                        if used[j]:
                            continue
                        diff = np.dot((coords[:, k] - crd[:, j]).T, (coords[:, k] - crd[:, j]))
                        if diff <= tol:
                            used[j] = 1
                            ifFound[k] = 1
                            break
                if np.sum(ifFound) == len(conf):
                    eFlag = 1
                    degens[i] += 1
                    break
            if eFlag:
                break
        if eFlag:
            continue
        etmp = []
        for i, W in enumerate(all_W):
            trans = np.zeros((4, 4))
            trans[0:3, 0:3] = W[:, :]
            trans[0:3, 3] = all_w[i][0:3]
            trans[3, :] = [0., 0., 0., 1.]
            new_coords = np.dot(trans, coords)
            for j in range(len(conf)):
                for k in range(3):
                    if np.isclose(new_coords[k, j], 0.):
                        new_coords[k, j] = 0.
                    elif np.isclose(new_coords[k, j], 1.):
                        new_coords[k, j] = 0.
                    elif new_coords[k, j] < 0:
                        new_coords[k, j] += 1.
                    elif new_coords[k, j] > 1:
                        new_coords[k, j] -= 1.
            etmp.append(np.unique(new_coords, axis=-1))
        all_eSites.append(etmp)
        unique.append(cind)
        degens.append(1)
    return unique, degens
