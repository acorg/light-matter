from unittest import TestCase

from light.performance import testArgs
from light.performance.affinity import AffinityMatrices, getScore
from light.performance.graphics import plot, plot3D
from light.performance.utils import makeOutputDir
from light.performance.data.polymerase import (
    Z_SCORES, BIT_SCORES, QUERIES, SUBJECTS, DATASET)
from light.performance.test.filesystem import FILESYSTEM_NAME


_AFFINITY = AffinityMatrices(QUERIES, SUBJECTS, returnDict=True,
                             returnAnalysis=True)


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

            allZScores = []
            allLmScores = []
            for queryId in Z_SCORES:
                zScores = []
                lmScores = []
                for subjectId in Z_SCORES[queryId]:
                    if queryId != subjectId:
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                        zScores.append(Z_SCORES[queryId][subjectId])
                allZScores.append(zScores)
                allLmScores.append(lmScores)
                plot(lmScores, zScores, queryId,
                     scoreTypeX, scoreTypeY, dirName)
            plot(allLmScores, allZScores, 'all',
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

            allBitScores = []
            allLmScores = []
            for queryId in BIT_SCORES:
                bitScores = []
                lmScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                allBitScores.append(bitScores)
                allLmScores.append(lmScores)
                plot(lmScores, bitScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)
            plot(allLmScores, allBitScores, 'all', scoreTypeX, scoreTypeY,
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

            allZScores = []
            allBitScores = []
            for queryId in BIT_SCORES:
                zScores = []
                bitScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                allZScores.append(zScores)
                allBitScores.append(bitScores)
                plot(bitScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)
            plot(allBitScores, allZScores, 'all', scoreTypeX, scoreTypeY,
                 dirName)


class TestBitScoreZScoreLightMatterScore3D(TestCase):

    def test3Dpolymerase(self):
        """
        Make a 3D plot of BLAST bit scores, Z scores, and light matter scores.
        """
        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(DATASET, parameterSet, '3d')
            allZScores = []
            allBitScores = []
            allLmScores = []
            for queryId in sorted(BIT_SCORES):
                zScores = []
                bitScores = []
                lmScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                allZScores.extend(zScores)
                allBitScores.extend(bitScores)
                allLmScores.extend(lmScores)

                plot3D(bitScores, zScores, lmScores, queryId,
                       'Bit score', 'Z score', 'Light matter score',
                       dirName, testArgs.interactive)

            # Plot all scores on one plot.
            plot3D(allBitScores, allZScores, allLmScores, 'all',
                   'Bit score', 'Z score', 'Light matter score',
                   dirName, testArgs.interactive)
