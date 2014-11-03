from unittest import TestCase

from light.statistics.bunyaviridae import Bunyaviridae, ALL_BUNYAVIRIDAE


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
        self.assertFalse(result)

    def testTwoPresent(self):
        """
        Must return True if two bunyavirus sequences are present.
        """
        sequence = 'TCATCACATGAATATATATATATATATCTCGTAGTATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertTrue(result)

    def testOnePresentBeginning(self):
        """
        Must return True if a bunyavirus sequence is present at the beginning.
        """
        sequence = 'TGTGTTTCATATATATATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertTrue(result)

    def testOnePresentEnd(self):
        """
        Must return True if a bunyavirus sequence is present at the end.
        """
        sequence = 'ATATATATATATATATGTGTTTC'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertTrue(result)

    def testOnePresentMiddle(self):
        """
        Must return True if a bunyavirus sequence is present in the middle.
        """
        sequence = 'ATATATATGTGTTTCTATATATA'
        statistic = Bunyaviridae()
        result = statistic.findBunyaviridae(sequence)
        self.assertTrue(result)

    def testAllDifferent(self):
        """
        All sequences must be distinct.
        """
        self.assertEqual(len(ALL_BUNYAVIRIDAE), len(set(ALL_BUNYAVIRIDAE)))
