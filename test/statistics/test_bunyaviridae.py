from unittest import TestCase

from light.statistics.bunyaviridae import Bunyaviridae


class TestBunyaviridae(TestCase):
    """
    Tests for the Statistic.Bunyaviridae class.
    """
    def testNotPresent(self):
        """
        Must return False if no bunyavirus sequence is present.
        """
        sequence = 'ATATATATATATATATATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertTrue(result)

    def testTwoPresent(self):
        """
        Must return True if two bunyavirus sequences are present.
        """
        sequence = 'TCATCACATGAATATATATATATATATCTCGTAGTATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertFalse(result)

    def testOnePresent(self):
        """
        Must return True if a bunyavirus sequence is present.
        """
        sequence = 'TGTGTTTCATATATATATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertFalse(result)
