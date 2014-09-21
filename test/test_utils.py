from unittest import TestCase

from light import utils
from dark.reads import Read, Reads


class TestStatistic(TestCase):
    """
    Tests for the Statistics class.
    """
    def testSupportedType(self):
        """
        Statistic.evaluate should raise an AttributeError, as SUPPORTED_TYPES
        is not defined
        """
        sequence = Read('id', 'ACGT')
        statistic = utils.Statistic()
        self.assertRaisesRegexp(Exception, 'Has no attribute SUPPORTED_TYPES',
                                statistic.evaluate, sequence)

    def testMinLength(self):
        """
        Statistic.evaluate should raise an AttributeError, as MIN_LENGTH
        is not defined
        """
        sequence = Read('id', 'ACGT')
        statistic = utils.Statistic()
        # adding a SUPPORTED_TYPES attribute, to allow test for MIN_LENGTH
        # to run.
        statistic.SUPPORTED_TYPES = True
        self.assertRaisesRegexp(Exception, 'Has no attribute MIN_LENGTH',
                                statistic.evaluate, sequence)


class TestHasAllBases(TestCase):
    """
    Tests for the HasAllBases subclass.
    """
    def testSupportedType(self):
        """
        If type of sequence not in SUPPORTED_TYPES, return False.
        """
        sequence = Read('id', 'ACGT')
        statistic = utils.HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertEqual(result, True)

    def testMinLength(self):
        """
        If the sequence is too short, return False.
        """
        sequence = Read('id', 'ACT')
        statistic = utils.HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertFalse(result)

    def testAllBases(self):
        """
        If the sequence has all possible bases, return True.
        """
        sequence = Read('id', 'ACTG')
        statistic = utils.HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertTrue(result)

    def testNotAllBases(self):
        """
        If the sequence does not contain all possible bases, return False.
        """
        sequence = Read('id', 'ACTT')
        statistic = utils.HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertFalse(result)

    def testHasName(self):
        """
        Test that the class has a name.
        """
        statistic = utils.HasAllBases()
        self.assertEqual('hasAllBases', statistic.NAME)


class TestRunStatistic(TestCase):
    """
    Tests for the runStastic function.
    """
    def testFunctionCall(self):
        statistic = utils.HasAllBases()
        read1 = Read('id1', 'ATGC')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        fastaReads = FastaReads()

        name, count = utils.runStatistic(statistic, fastaReads)
        self.assertEqual('hasAllBases', name)
        self.assertEqual(1, count)
