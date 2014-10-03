from unittest import TestCase

from light.statistics.alpha_helix import AlphaHelix
from dark.reads import Read


class TestAlphaHelix(TestCase):
    """
    Tests for the Statistic.AlphaHelix class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        sequence = 'ASDGEBHSDTDSCV'
        result = AlphaHelix.convertAAToHydrophobicHydrophilic(sequence)
        self.assertEqual('IOOIOIOOIOOII', result)

    def test_evaluateNoHelix(self):
        """
        The AlphaHelix._evaluate method should return False when no alpha helix
        is present.
        """
        sequence = Read('id', 'ASDGEBHSDTDSCV')
        statistic = AlphaHelix()
        result = statistic.evaluate(sequence)
        self.assertFalse(result)

    def test_evaluateHelix(self):
        """
        The AlphaHelix._evaluate method should return True when an alpha helix
        is present.
        """
        sequence = Read('id', 'FFFRFFFRFFFRFFFR')
        statistic = AlphaHelix()
        result = statistic.evaluate(sequence)
        self.assertTrue(result)
