from unittest import TestCase

from .filesystem import FILESYSTEM_NAME

from light.performance import testArgs
from light.performance.affinity import AffinityMatrices
from light.performance.graphics import plot, plot3D
from light.performance.utils import makeOutputDir
from light.performance.data._4mtp_a import (
    Z_SCORES, BIT_SCORES, QUERIES, SUBJECTS, DATASET)

_AFFINITY = AffinityMatrices(QUERIES, SUBJECTS, symmetric=False,
                             returnDict=True)


class TestLightMatterZScoreCorrelation(TestCase):

    def testPlotsMTP(self):
        """
        Examine the correlation between our scores and Z scores.
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
                lmScores = []
                zScores = []
                for subjectId in Z_SCORES[queryId]:
                    if queryId != subjectId:
                        zScores.append(Z_SCORES[queryId])
                        lmScores.append(affinity[queryId][subjectId])
                plot(lmScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestLightMatterBitScoreCorrelation(TestCase):

    def testPlotsMTP(self):
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
                lmScores = []
                bitScores = []
                for subjectId in BIT_SCORES[queryId]:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId])
                        lmScores.append(affinity[queryId][subjectId])
                plot(lmScores, bitScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreCorrelation(TestCase):

    def testPlotsMTP(self):
        """
        Examine the correlation between Z scores and blast bit scores.
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
                        bitScores.append(BIT_SCORES[queryId])
                        zScores.append(Z_SCORES[queryId])
                plot(bitScores, zScores, queryId, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreLightMatterScore3D(TestCase):

    def test2HLA(self):
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
