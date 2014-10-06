from unittest import TestCase

from light.statistics.all_bases import HasAllBases
from dark.reads import Read


class TestHasAllBases(TestCase):
    """
    Tests for the HasAllBases subclass.
    """
    def testSupportedType(self):
        """
        If type of sequence not in SUPPORTED_TYPES, return False.
        """
        sequence = Read('id', 'ACGT')
        statistic = HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertTrue(result)

    def testMinLength(self):
        """
        If the sequence is too short, return False.
        """
        sequence = Read('id', 'ACT')
        statistic = HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertFalse(result)

    def testAllBases(self):
        """
        If the sequence has all possible bases, return True.
        """
        sequence = Read('id', 'ACTG')
        statistic = HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertTrue(result)

    def testNotAllBases(self):
        """
        If the sequence does not contain all possible bases, return False.
        """
        sequence = Read('id', 'ACTT')
        statistic = HasAllBases()
        result = statistic.evaluate(sequence)
        self.assertFalse(result)

    def testHasName(self):
        """
        Test that the class has a name.
        """
        statistic = HasAllBases()
        self.assertEqual('hasAllBases', statistic.NAME)
