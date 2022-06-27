import numpy as np


def stepped_enumerate(xs, start=0, step=1):
    """
    Adaptation to the builtin `enumerate` with a `step` arg.
    """
    for x in xs:
        yield (start, x)
        start += step


def get_complement(mylist, idx):
    """
    Returns the sublist of `mylist` consisting of the complement of elements indexed by `idx`.
    """
    myarray = np.array(mylist)
    mask = np.full(len(mylist), False)
    mask[idx] = True
    complement = list(myarray[~mask])

    return complement
