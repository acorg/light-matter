from unittest import TestCase

from light.performance import testArgs
from light.performance.affinity import AffinityMatrices
from light.performance.utils import (
    plot, plot3D, makeOutputDir, FILESYSTEM_NAME)
from light.performance.data.polymerase import (
    Z_SCORES, BIT_SCORES, QUERIES, SUBJECTS, DATASET)


_AFFINITY = AffinityMatrices(QUERIES, SUBJECTS)


class TestLightMatterZScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between light matter scores and Z scores.
        """
        scoreTypeX = 'Light matter score'
        scoreTypeY = 'Z score'

        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(
                DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))
            for queryId in Z_SCORES:
                zScores = []
                lmScores = []
                for subjectId in Z_SCORES[queryId]:
                    if queryId != subjectId:
                        lmScores.append(affinity[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                plot(lmScores, zScores, queryId,
                     scoreTypeX, scoreTypeY, dirName)


class TestLightMatterBitScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between our scores and blast bit scores.
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
            for queryId in BIT_SCORES:
                zScores = []
                lmScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        lmScores.append(affinity[queryId][subjectId])
                        zScores.append(BIT_SCORES[queryId][subjectId])
                plot(lmScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between BLAST bit scores and Z scores.
        """
        scoreTypeX = 'Bit score'
        scoreTypeY = 'Z score'

        for parameterSet in testArgs.parameterSets:
            dirName = makeOutputDir(
                DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))

            for queryId in BIT_SCORES:
                zScores = []
                bitScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                plot(bitScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreLightMatterScore3D(TestCase):

    def test3Dpolymerase(self):
        """
        Make a 3D plot of BLAST bit scores, Z scores, and light matter scores.
        """
        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(DATASET, parameterSet, '3d')
            for queryId in sorted(BIT_SCORES):
                zScores = []
                bitScores = []
                lmScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                        lmScores.append(affinity[queryId][subjectId])
                plot3D(bitScores, zScores, lmScores, queryId,
                       'Bit score', 'Z score', 'Light matter score',
                       dirName, testArgs.interactive)
