from unittest import TestCase

from light.statistics import Statistic, runStatistic
from light.statistics.all_bases import HasAllBases
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
    def testFunctionCall(self):
        statistic = HasAllBases()
        read1 = Read('id1', 'ATGC')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        fastaReads = FastaReads()

        name, count = runStatistic(statistic, fastaReads)
        self.assertEqual('hasAllBases', name)
        self.assertEqual(1, count)
