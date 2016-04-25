from os.path import dirname, join
from unittest import TestCase

from dark.reads import SSAARead, SSAAReadWithX, Reads
from dark.fasta_ss import SSFastaReads

import light
from light.database import Database
from light.performance import testArgs
from light.performance.affinity import affinityMatrix
from light.performance.data.hla import BIT_SCORES, Z_SCORES
from light.performance.test.utils import plot


# Create a singleton affinity matrix of lm scores for all sequences.
_SUBJECT = Reads()
_SUBJECT.add(SSAARead('2HLA:A:sequence',
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

LM_SCORES = affinityMatrix(_QUERIES, subjects=_SUBJECT, symmetric=False,
                           database=Database(testArgs.dbParams),
                           findParams=testArgs.findParams, returnDict=True)


class TestZScoreCorrelation(TestCase):

    def testPlotsHLA(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        lmScores = []
        zScores = []
        for query in _QUERIES:
            assert query.id.endswith(':sequence')
            id_ = query.id[:-9]
            zScores.append(Z_SCORES[id_])
            lmScores.append(LM_SCORES[query.id]['2HLA:A:sequence'])
        plot(lmScores, zScores, '2HLA:A', 'Light matter score', 'Z score',
             'hla')


class TestBitScoreCorrelation(TestCase):

    def testPlotsHLA(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        lmScores = []
        bitScores = []
        for query in _QUERIES:
            assert query.id.endswith(':sequence')
            id_ = query.id[:-9]
            bitScores.append(BIT_SCORES[id_])
            lmScores.append(LM_SCORES[query.id]['2HLA:A:sequence'])
        plot(lmScores, bitScores, '2HLA:A', 'Light matter score',
             'Bit score', 'hla')


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlotsHLA(self):
        """
        Examine the correlation between Z scores and blast bit scores.
        """
        zScores = []
        bitScores = []
        for query in _QUERIES:
            assert query.id.endswith(':sequence')
            id_ = query.id[:-9]
            bitScores.append(BIT_SCORES[id_])
            zScores.append(Z_SCORES[id_])
        plot(bitScores, zScores, '2HLA:A', 'Bit score', 'Z score', 'hla')
