from unittest import TestCase

from light.performance import testArgs
from light.performance.affinity import AffinityMatrices, getScore
from light.performance.data.pdb_2hla_a_against_polymerase import (
    BIT_SCORES, QUERIES, SUBJECTS, DATASET)
from light.performance.graphics import plot
from light.performance.test.filesystem import FILESYSTEM_NAME
from light.performance.utils import makeOutputDir

_AFFINITY = AffinityMatrices(QUERIES, SUBJECTS, symmetric=False,
                             computeDiagonal=True, returnDict=True,
                             returnAnalysis=True)


class TestLightMatterBitScoreCorrelation(TestCase):

    def testNonMatchingPlots2HLAAgainstPolymerase(self):
        """
        Examine the correlation between light matter scores and BLAST bit
        scores.
        """
        scoreTypeX = 'Light matter score'
        scoreTypeY = 'Bit score'

        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(
                DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))
            for subject in SUBJECTS:
                lmScores = []
                bitScores = []
                for query in QUERIES:
                    if query.id != subject.id:
                        bitScores.append(BIT_SCORES[query.id][subject.id])
                        lmScore = getScore(affinity, query.id, subject.id)
                        lmScores.append(lmScore)
                plot(lmScores, bitScores, subject.id, scoreTypeX, scoreTypeY,
                     dirName)
