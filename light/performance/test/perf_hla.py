from os.path import dirname, join
from unittest import TestCase

from dark.reads import SSAARead, SSAAReadWithX, Reads
from dark.fasta_ss import SSFastaReads

import light
from light.performance import testArgs
from light.performance.affinity import AffinityMatrices
from light.performance.data.hla import BIT_SCORES, Z_SCORES
from light.performance.utils import (
    plot, plot3D, makeOutputDir, FILESYSTEM_NAME)

_DATASET = 'hla'
_SUBJECT_NAME = '2HLA:A'

_SUBJECT = Reads()
_SUBJECT.add(SSAARead(_SUBJECT_NAME,
                      'GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQE'
                      'GPEYWDRNTRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQMMYGCDVGSDGRFL'
                      'RGYRQDAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQWRAYLEGTCV'
                      'EWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRD'
                      'GEDQTQDTELVETRPAGDGTFQKWVAVVVPSGQEQRYTCHVQHEGLPKPL',
                      '--EEEEEEEEEE--TTSS--EEEEEEEETTEEEEEEETTSTT-S-EE-SHHHHTS'
                      '-HHHHHHHHHHHHHHHHHHHHHHHHHHHHTT--TTS--EEEEEEEEEE-TTS-EE'
                      'EEEEEEEETTEEEEEE-TTSS-EEESSHHHHHHHHHHHHTTTHHHHHHHHHTHHH'
                      'HHHHHHHHHHHHHHT--B--EEEEEEEE-SSSEEEEEEEEEEEBSS-EEEEEEET'
                      'TEEE-TTEEE---EE-SSS-EEEEEEEEEETT-GGGEEEEEEETTB-S--'))

_QUERIES = list(SSFastaReads(
    join(dirname(light.__file__),
         'performance', 'data', 'hla-queries-structure.fasta'),
    readClass=SSAAReadWithX))

_AFFINITY = AffinityMatrices(_QUERIES, _SUBJECT, symmetric=False)


class TestZScoreCorrelation(TestCase):

    def testPlotsHLA(self):
        """
        Examine the correlation between light matter scores and Z scores.
        """
        scoreTypeX = 'Light matter score'
        scoreTypeY = 'Z score'

        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(
                _DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))
            lmScores = []
            zScores = []
            for query in _QUERIES:
                zScores.append(Z_SCORES[query.id])
                lmScores.append(affinity[query.id][_SUBJECT_NAME])
            plot(lmScores, zScores, _SUBJECT_NAME, scoreTypeX, scoreTypeY,
                 dirName)


class TestBitScoreCorrelation(TestCase):

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
                _DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))
            lmScores = []
            bitScores = []
            for query in _QUERIES:
                bitScores.append(BIT_SCORES[query.id])
                lmScores.append(affinity[query.id][_SUBJECT_NAME])
            plot(lmScores, bitScores, _SUBJECT_NAME, scoreTypeX,
                 scoreTypeY, dirName)


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlotsHLA(self):
        """
        Examine the correlation between BLAST bit scores and Z scores.
        """
        scoreTypeX = 'Bit score'
        scoreTypeY = 'Z score'

        for parameterSet in testArgs.parameterSets:
            dirName = makeOutputDir(
                _DATASET,
                parameterSet,
                '%s-%s' % (FILESYSTEM_NAME[scoreTypeX],
                           FILESYSTEM_NAME[scoreTypeY]))
            zScores = []
            bitScores = []
            for query in _QUERIES:
                bitScores.append(BIT_SCORES[query.id])
                zScores.append(Z_SCORES[query.id])
            plot(bitScores, zScores, _SUBJECT_NAME, scoreTypeX, scoreTypeY,
                 dirName)


class TestZScoreBitScoreLMScore3D(TestCase):

    def test3DHLA(self):
        """
        Make a 3D plot of BLAST bit scores, Z scores, and light matter scores.
        """
        for parameterSet in testArgs.parameterSets:
            affinity = _AFFINITY[parameterSet]
            dirName = makeOutputDir(_DATASET, parameterSet, '3d')
            for queryId in sorted(BIT_SCORES):
                zScores = []
                bitScores = []
                lmScores = []
                for subjectId in BIT_SCORES:
                    if queryId != subjectId:
                        bitScores.append(BIT_SCORES[queryId][subjectId])
                        zScores.append(Z_SCORES[queryId][subjectId])
                        lmScores.append(affinity[queryId][subjectId])
                plot3D(bitScores, zScores, lmScores, queryId,
                       'Bit score', 'Z score', 'Light matter score',
                       dirName, testArgs.interactive)
