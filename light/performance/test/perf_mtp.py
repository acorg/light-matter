from os.path import dirname, join
from unittest import TestCase

from dark.reads import SSAARead, SSAAReadWithX, Reads
from dark.fasta_ss import SSFastaReads

import light
from light.database import Database
from light.performance import testArgs
from light.performance.affinity import affinityMatrix
from light.performance.data.mtp import BIT_SCORES, Z_SCORES
from light.performance.test.utils import plot

# Create a singleton affinity matrix of lm scores for all sequences.
_SUBJECT = Reads()
_SUBJECT.add(SSAARead('4MTP:A:sequence',
                      'VHSNQEKIKKRIQKLKEEFATTWHKDPEHPYRTWTYHGSYEVKATGSASSLVNGV'
                      'VKLMSKPWDAIANVTTMAMTDTTPFGQQRVFKEKVDTKAPEPPAGVREVLNETTN'
                      'WLWAHLSREKRPRLCTKEEFIKKVNSNAALGAVFAEQNQWSTAREAVNDPRFWEM'
                      'VDEERENHLRGECHTCIYNMMGKREKKPGEFGKAKGSRAIWFMWLGARYLEFEAL'
                      'GFLNEDHWLSRENSGGGVEGSGVQKLGYILRDIAGKQGGKMYADDTAGWDTRITR'
                      'TDLENEAKVLELLDGEHRMLARAIIELTYRHKVVKVMRPAAEGKTVMDVISREDQ'
                      'RGSGQVVTYALNTFTNIAVQLVRLMEAEGVIGPQHLEQLPRKNKIAVRTWLFENG'
                      'EERVTRMAISGDDCVVKPLDDRFATALHFLNAMSKVRKDIQEWKPSHGWHDWQQV'
                      'PFCSNHFQEIVMKDGRSIVVPCRGQDELIGRARISPGAGWNVKDTACLAKAYAQM'
                      'WLLLYFHRRDLRLMANAICSAVPVDWVPTGRTSWSIHSKGEWMTTEDMLQVWNRV'
                      'WIEENEWMMDKTPIASWTDVPYVGKREDIWCGSLIGTRSRATWAENIYAAINQVR'
                      'AVIGKENYVDYMTSLRRYEDVLIQEDRVI',
                      '------TTHHHHHHHHHHTTTT-B---S---SSSEEEEEEE----------B-HH'
                      'HHHH-GGGGS-HHHHT---EE-SHHHHHHHHHHHTS-------HHHHHHHHHHHH'
                      'HHHHHHTTT-------HHHHHHHH-------------------HHHHT-TTHHHH'
                      'HHHHHHHHHHT--SS--EEEEEEE----SSSS-----EEEEE--HHHHHHHHHHH'
                      'THHHHTTTTSHHHHSSB-TT--HHHHHHHHHHHHHSSSS-EE---BTTGGGG--H'
                      'HHHHHHGGGGGG--THHHHHHHHHHHHTTTSEEEEEEEEETTTEEEEEEEEESS-'
                      '--S-STTHHHHHHHHHHHHHHHHHHHHHTSS-GGGTTS--HHHHHHHHHHHHHHH'
                      'HHHGGGEEEETTEEEE--SSGGGGG--HHHHHTT--BSSS-TTS---EES-GGG-'
                      '-BTTBEEEEEE-TTS-EEEEEE--HHHHHHHHHB-------HHHHHHHHHHHHHH'
                      'HHHHSTTSHHHHHHHHHHHHTS-TT----S-S---TT---TTSSSS-HHHHHHHH'
                      'HTTS-TT--------SGGGS----HHHHHHTT--TT-HHHHHHHHTHHHHHHHHH'
                      'HHH-S------------------------'))

_QUERIES = list(SSFastaReads(
    join(dirname(light.__file__),
         'performance', 'data', 'mtp-queries-structure.fasta'),
    readClass=SSAAReadWithX))

LM_SCORES = affinityMatrix(_QUERIES, subjects=_SUBJECT, symmetric=False,
                           database=Database(testArgs.dbParams),
                           findParams=testArgs.findParams, returnDict=True)


class TestZScoreCorrelation(TestCase):

    def testPlotsMTP(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        lmScores = []
        zScores = []
        for query in _QUERIES:
            assert query.id.endswith(':sequence')
            id_ = query.id[:-9]
            zScores.append(Z_SCORES[id_])
            lmScores.append(LM_SCORES[query.id]['4MTP:A:sequence'])
        plot(lmScores, zScores, '4MPT:A', 'Light matter score', 'Z score',
             'mtp')


class TestBitScoreCorrelation(TestCase):

    def testPlotsMTP(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        lmScores = []
        bitScores = []
        for query in _QUERIES:
            assert query.id.endswith(':sequence')
            id_ = query.id[:-9]
            bitScores.append(BIT_SCORES[id_])
            lmScores.append(LM_SCORES[query.id]['4MTP:A:sequence'])
        plot(lmScores, bitScores, '4MPT:A', 'Light matter score',
             'Bit score', 'mtp')


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlotsMTP(self):
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
        plot(bitScores, zScores, '4MPT:A', 'Bit score', 'Z score', 'mtp')
