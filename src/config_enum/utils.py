"""General helper functions of config_enum (names as in the original script)."""
import operator as op
import random
import re
import string
from functools import reduce
from subprocess import PIPE, Popen

import numpy as np


def Popen4(cmd):
    """Run a shell command; returns (list of non-empty output lines, stderr)."""
    proc = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out, err = proc.communicate()
    out2 = [x for x in out.decode(errors="replace").split("\n") if x != ""]
    return out2, err


def hashid(n=6):
    return ''.join([random.choice(string.ascii_lowercase + string.digits) for _ in range(n)])


def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect."""
    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def cartesian(arrays, out=None):  # Not required anymore.
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))[:3]
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6]])
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = int(np.prod([x.size for x in arrays]))
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n // arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
    return out


def combinations2(iterable, r):  # from itertools
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(list(range(r))):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)


def ncr(n, r):  # gives combination.
    r = min(r, n - r)
    if r < 0:
        return 0
    if r == 0:
        return 1
    numer = reduce(op.mul, list(range(n, n - r, -1)))
    denom = reduce(op.mul, list(range(1, r + 1)))
    return numer // denom


def sort_dict_by_val(data, index=0):
    """Keys of a dictionary sorted by their values (ascending).

    index: index to use if values of the dictionary is not scalar.
    """
    keys = list(data.keys())
    if len(keys) <= 1:
        return keys

    skeys = [keys[0]]  # sorted keys
    for i in range(len(keys)):
        key = keys[i]
        val = data[key]
        try:
            val = val[index]
        except (TypeError, IndexError):
            pass
        flag = 0
        for j in range(len(skeys)):
            key2 = skeys[j]
            val2 = data[key2]
            try:
                val2 = val2[index]
            except (TypeError, IndexError):
                pass
            if val <= val2 and i != 0:
                skeys.insert(j, key)
                flag = 1
                break
        if not flag and key not in skeys:
            skeys.append(key)  # not to miss the largest element.

    return skeys


def parse_ranges(entries, sep=":"):
    """'0:4 24' (or ['0:4', '24']) -> [0, 1, 2, 3, 4, 24] (ranges include both ends)."""
    if isinstance(entries, str):
        entries = entries.replace(",", " ").split()
    vals = []
    for y in entries:
        for z in str(y).replace(",", " ").split():
            x = z.split(sep)
            if len(x) > 1:
                vals.extend(list(range(int(x[0]), int(x[1]) + 1)))
            else:
                vals.append(int(x[0]))
    return vals
