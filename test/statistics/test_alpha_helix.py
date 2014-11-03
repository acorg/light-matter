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
        self.assertEqual('OIIOIOOIIOIIOO', result)

    def testFindWithHelix(self):
        """
        The find method must return the right position of the helix.
        """
        sequence = 'OIIIOIIIOIOIOIOIOIOIOIOIOIOIOIOOIIIOIIIO'
        statistic = AlphaHelix()
        result = statistic.find(sequence)
        self.assertEqual([0, 31], result)

    def testFindWithoutHelix(self):
        """
        The find method must return false when no helix is there.
        """
        sequence = 'OIOIOIOIOIOIOIOIOIOO'
        statistic = AlphaHelix()
        result = statistic.find(sequence)
        self.assertEqual([], result)

    def test_evaluateNoHelix(self):
        """
        The AlphaHelix._evaluate method should return False when no alpha helix
        is present.
        """
        sequence = Read('id', 'ASDGEAHSDTDSCV')
        statistic = AlphaHelix()
        result = statistic.evaluate(sequence)
        self.assertEqual(False, result)

    def test_evaluateHelix(self):
        """
        The AlphaHelix._evaluate method should return True when an alpha helix
        is present.
        """
        sequence = Read('id', 'RRRFRRRFRRRFRRRFFFFFFF', type='aa')
        statistic = AlphaHelix()
        result = statistic._evaluate(sequence)
        self.assertEqual(1, result)

    def test_evaluateHelixDistanceThreeHelices(self):
        """
        The calculateDistance function needs to return the right distance
        between alpha helices.
        """
        sequence = Read('id', 'RRRFRRRFRRRFRRRFFFFFFFRRRFRRRFRRRFRRR',
                        type='aa')
        statistic = AlphaHelix()
        result = statistic._evaluate(sequence, distances=True)
        self.assertEqual([4, 14, 4], result)

    def test_evaluateHelixDistanceOneHelix(self):
        """
        The calculateDistance function needs to return the right distance
        between alpha helices.
        """
        sequence = Read('id', 'RRRFRRRFRRRF',
                        type='aa')
        statistic = AlphaHelix()
        result = statistic._evaluate(sequence, distances=True)
        self.assertFalse(result)
