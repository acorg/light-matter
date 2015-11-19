from itertools import permutations
import numpy as np

from light.performance import affinity


def ultrametric(sequenceFileOrMatrix, findParams=None, **kwargs):
    """
    Test whether ultrametricity is satisfied for a distance matrix.
    Ultrametricity is satisfied when: d(A, C) <= max(d(A, B), d(B, C)) for any
    three scores (from sequence comparisons) A, B, C.

    @param sequenceFileOrMatrix: Either a C{str} file name of a file
        containing sequences or a distance matrix as returned from
        C{light.performance.affinity}.
    @param findParams: A C{light.parameters.FindParameters} instance.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.

    @return: A generator which returns non-ultrametric triplets.
    """
    if isinstance(sequenceFileOrMatrix, np.ndarray):
        matrix = sequenceFileOrMatrix

    else:
        matrix = affinity.affinityMatrix(sequenceFileOrMatrix, findParams,
                                         **kwargs)

    for a, b, c in permutations(range(len(matrix)), 3):
        if matrix[a][c] < max(matrix[a][b], matrix[b][c]):
            yield a, b, c
