from os.path import dirname, join
from unittest import TestCase

from dark.fasta_ss import SSFastaReads

import light
from light.database import Database
from light.performance import testArgs
from light.performance.affinity import affinityMatrix
from light.performance.polymerase import Z_SCORES, BIT_SCORES
from light.performance.test.utils import plot

# Create a singleton affinity matrix of lm scores for all polymerase
# sequences.
_QUERIES = list(SSFastaReads(
    join(dirname(light.__file__),
         'performance', 'data', 'polymerase.fasta')))
_AFFINITY = affinityMatrix(_QUERIES, database=Database(testArgs.dbParams),
                           findParams=testArgs.findParams, returnDict=True)


class TestZScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        for queryId in Z_SCORES:
            zScores = []
            lmScores = []
            for subjectId in Z_SCORES:
                if queryId != subjectId:
                    lmScores.append(_AFFINITY[queryId][subjectId])
                    zScores.append(Z_SCORES[queryId][subjectId])
            plot(lmScores, zScores, queryId, 'Light matter score', 'Z score',
                 'polymerase')


class TestBitScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        for queryId in BIT_SCORES:
            zScores = []
            lmScores = []
            for subjectId in BIT_SCORES:
                if queryId != subjectId:
                    lmScores.append(_AFFINITY[queryId][subjectId])
                    zScores.append(BIT_SCORES[queryId][subjectId])
            plot(lmScores, zScores, queryId, 'Light matter score', 'Bit score',
                 'polymerase')


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between Z scores and blast bit scores.
        """
        for queryId in BIT_SCORES:
            zScores = []
            bitScores = []
            for subjectId in BIT_SCORES:
                if queryId != subjectId:
                    bitScores.append(BIT_SCORES[queryId][subjectId])
                    zScores.append(Z_SCORES[queryId][subjectId])
            plot(bitScores, zScores, queryId, 'Bit score', 'Z score',
                 'polymerase')
