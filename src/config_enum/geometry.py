"""Distances and cell helpers (names as in the original script)."""
import math

import numpy as np


def dist(arr1, arr2):
    # Distance btw two atoms. Works with two 1x3 arrays.
    if len(arr1) == 4:
        arr1 = arr1[1:4]
    if len(arr2) == 4:
        arr2 = arr2[1:4]
    return math.sqrt((arr1[0] - arr2[0])**2 + (arr1[1] - arr2[1])**2 + (arr1[2] - arr2[2])**2)


def minDist(arr1, arr2, latt):
    # Finds the minimum distance from the images applying the PBC. (Minimum image convention).
    if len(arr1) == 4:
        arr1 = arr1[1:4]
    if len(arr2) == 4:
        arr2 = arr2[1:4]
    minD = dist(arr1, arr2)
    for i in range(-1, 2, 1):
        for j in range(-1, 2, 1):
            for k in range(-1, 2, 1):
                newarr2 = [0.0, 0.0, 0.0]
                newarr2[0] = arr2[0] + i * latt[0][0] + j * latt[1][0] + k * latt[2][0]
                newarr2[1] = arr2[1] + i * latt[0][1] + j * latt[1][1] + k * latt[2][1]
                newarr2[2] = arr2[2] + i * latt[0][2] + j * latt[1][2] + k * latt[2][2]
                currD = dist(arr1, newarr2)
                if currD < minD:
                    minD = currD
    return minD


def volume(cell):
    return np.abs(np.dot(cell[2], np.cross(cell[0], cell[1])))


def lenCell(cell):
    return [np.linalg.norm(v) for v in cell]


def lenAngleCell(cell, radians=False):  # taken from ase.geometry.cell_to_cellpar
    lengths = [np.linalg.norm(v) for v in cell]
    angles = []
    for i in range(3):
        j = i - 1
        k = i - 2
        ll = lengths[j] * lengths[k]
        if ll > 1e-16:
            x = np.dot(cell[j], cell[k]) / ll
            angle = 180.0 / np.pi * np.arccos(x)
        else:
            angle = 90.0
        angles.append(angle)
    if radians:
        angles = [angle * np.pi / 180 for angle in angles]
    return np.array(lengths + angles)


def compCells(cell1, cell2):  # Compare cells.
    len1 = lenCell(cell1)
    len2 = lenCell(cell2)
    return tuple([int(round(len1[i] / len2[i])) for i in range(3)])
