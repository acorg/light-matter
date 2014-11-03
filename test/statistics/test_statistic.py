from unittest import TestCase

from light.statistics import Statistic, runStatistic, find
from light.statistics.all_bases import HasAllBases
from light.statistics.bunyaviridae import Bunyaviridae
from light.statistics.alpha_helix import AlphaHelix
from dark.reads import Read, Reads


class TestStatistic(TestCase):
    """
    Tests for the light.statistics.Statistic class.
    """
    def testSupportedType(self):
        """
        Statistic.evaluate should raise an AttributeError, as SUPPORTED_TYPES
        is not defined
        """
        sequence = Read('id', 'ACGT')
        statistic = Statistic()
        self.assertRaisesRegexp(Exception, 'Has no attribute SUPPORTED_TYPES',
                                statistic.evaluate, sequence)

    def testMinLength(self):
        """
        Statistic.evaluate should raise an AttributeError, as MIN_LENGTH
        is not defined
        """
        sequence = Read('id', 'ACGT')
        statistic = Statistic()
        # adding a SUPPORTED_TYPES attribute, to allow test for MIN_LENGTH
        # to run.
        statistic.SUPPORTED_TYPES = True
        self.assertRaisesRegexp(Exception, 'Has no attribute MIN_LENGTH',
                                statistic.evaluate, sequence)


class TestRunStatistic(TestCase):
    """
    Tests for the light.statistics.runStatistic function.
    """
    def testRunStatisticFunctionCall(self):
        statistic = HasAllBases()
        read1 = Read('id1', 'ATGC')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        fastaReads = FastaReads()

        count = runStatistic(statistic, fastaReads)
        self.assertEqual(1, count)


class TestFind(TestCase):
    """
    Tests for the light.statistics.find function.
    """
    def testFindAll(self):
        """
        Function should return all statistics, if asked for all statistics.
        """
        result = find()
        self.assertEqual([HasAllBases, AlphaHelix, Bunyaviridae], list(result))

    def testFindDNA(self):
        """
        Function should return statistics for dna, if asked for dna statistics.
        """
        result = find(worksOn='dna')
        self.assertEqual([HasAllBases, Bunyaviridae], list(result))

    def testFindProtein(self):
        """
        Function should return statistics for proteins, if asked for protein
        statistics.
        """
        result = find(worksOn='protein')
        self.assertEqual([], list(result))
