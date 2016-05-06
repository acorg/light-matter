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
        scoreTypeX = 'Secondary structure as measured by\nLight matter score'
        scoreTypeY = 'Tertiary structure as measured by\nZ score'

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
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                        zScores.append(Z_SCORES[queryId][subjectId])
                plot(lmScores, zScores, queryId,
                     scoreTypeX, scoreTypeY, dirName)


class TestLightMatterBitScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        scoreTypeX = 'Secondary structure as measured by\nLight matter score'
        scoreTypeY = 'Primary structure as measured by\nBit score'

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
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                        zScores.append(BIT_SCORES[queryId][subjectId])
                plot(lmScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreCorrelation(TestCase):

    def testPlotsPolymerase(self):
        """
        Examine the correlation between BLAST bit scores and Z scores.
        """
        scoreTypeX = 'Primary structure as measured by\nBit score'
        scoreTypeY = 'Tertiary structure as measured by\nZ score'

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
        scoreTypeX = 'Primary structure as measured by\nBit score'
        scoreTypeY = 'Tertiary structure as measured by\nZ score'
        scoreTypeZ = 'Secondary structure as measured by\nlight matter score'

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
                        lmScore = getScore(affinity, queryId, subjectId)
                        lmScores.append(lmScore)
                plot3D(bitScores, zScores, lmScores, queryId,
                       scoreTypeX, scoreTypeY, scoreTypeZ, dirName,
                       testArgs.interactive)
