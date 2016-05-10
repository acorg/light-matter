from unittest import TestCase

from light.performance import testArgs
from light.performance.affinity import AffinityMatrices, getScore
from light.performance.graphics import plot, plot3D
from light.performance.utils import makeOutputDir
from light.performance.data.pdb_2hla_a import (
    Z_SCORES, BIT_SCORES, QUERIES, SUBJECTS, DATASET)
from light.performance.test.filesystem import FILESYSTEM_NAME

_AFFINITY = AffinityMatrices(QUERIES, SUBJECTS, symmetric=False,
                             computeDiagonal=True, returnDict=True,
                             returnAnalysis=True)


class TestLightMatterZScoreCorrelation(TestCase):

    def testPlotsHLA(self):
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
            for subject in SUBJECTS:
                lmScores = []
                zScores = []
                for query in QUERIES:
                    if query.id != subject.id:
                        zScores.append(Z_SCORES[query.id][subject.id])
                        lmScore = getScore(affinity, query.id, subject.id)
                        lmScores.append(lmScore)
                plot(lmScores, zScores, subject.id, scoreTypeX, scoreTypeY,
                     dirName)


class TestLightMatterBitScoreCorrelation(TestCase):

    def testPlotsHLA(self):
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


class TestBitScoreZScoreCorrelation(TestCase):

    def testPlotsHLA(self):
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
            for subject in SUBJECTS:
                bitScores = []
                zScores = []
                for query in QUERIES:
                    if query.id != subject.id:
                        zScores.append(Z_SCORES[query.id][subject.id])
                        bitScores.append(BIT_SCORES[query.id][subject.id])
                plot(bitScores, zScores, subject.id, scoreTypeX, scoreTypeY,
                     dirName)


class TestBitScoreZScoreLightMatterScore3D(TestCase):

    def testPlotsHLA3D(self):
        """
        Make a 3D plot of BLAST bit scores, Z scores, and light matter scores.
        """
        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(DATASET, parameterSet, '3d')
            for subject in SUBJECTS:
                zScores = []
                bitScores = []
                lmScores = []
                for query in QUERIES:
                    if query.id != subject.id:
                        bitScores.append(BIT_SCORES[query.id][subject.id])
                        zScores.append(Z_SCORES[query.id][subject.id])
                        lmScore = getScore(affinity, query.id, subject.id)
                        lmScores.append(lmScore)
                plot3D(bitScores, zScores, lmScores, subject.id,
                       'Bit score', 'Z score', 'Light matter score',
                       dirName, testArgs.interactive)
