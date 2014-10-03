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
        sequence = 'ASDGEAHSDTDSCV'
        statistic = AlphaHelix()
        result = statistic.convertAAToHydrophobicHydrophilic(sequence)
        self.assertEqual(result, 'OIIOIOOIIOIIOO')

    def testFindWithHelix(self):
        """
        The find method must return the right position of the helix.
        """
        sequence = 'OIIIOIIIOIOIOIOIOIOIOIOIOIOIOIOOIIIOIIIO'
        statistic = AlphaHelix()
        result = statistic.find(sequence)
        self.assertEqual(result, [0, 31])

    def testFindWithoutHelix(self):
        """
        The find method must return false when no helix is there.
        """
        pass

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
        sequence = Read('id', 'RRRFRRRFRRRFRRRFFFFFFF')
        statistic = AlphaHelix()
        result = statistic.evaluate(sequence)
        self.assertTrue(result)
