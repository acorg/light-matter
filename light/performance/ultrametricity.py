from itertools import combinations
import numpy as np

from light.performance import affinity


def ultrametric(sequenceFileOrMatrix, findParams=None, **kwargs):
    """
    Test whether ultrametricity is satisfied for a distance matrix.
    Ultrametricity is satisfied when: dAC <= max(dAB, dBC)

    @param sequenceFileOrMatrix: Either a C{str} file name of a file
        containing sequences or a distance matrix as returned from
        C{light.performance.affinity}.
    @param findParams: A C{light.parameters.FindParameters} instance.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    if isinstance(sequenceFileOrMatrix, np.ndarray):
        matrix = sequenceFileOrMatrix

    else:
        matrix = affinity.affinityMatrix(sequenceFileOrMatrix, findParams,
                                         **kwargs)

    # get all possible triplets
    allTriplets = list(combinations(range(len(matrix)), 3))

    ultrametric = []

    for triplet in allTriplets:
        # get scores
        a = triplet[0]
        b = triplet[1]
        c = triplet[2]

        dac = matrix[a][c]
        dab = matrix[a][b]
        dbc = matrix[b][c]
        score = (dac <= max(dab, dbc))
        ultrametric.append(score)

    true = ultrametric.count(True)

    print('%d of %d triplets satisfy the ultrametricity criterion.' % (true,
          len(ultrametric)))
    return ultrametric, allTriplets
