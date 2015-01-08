from unittest import TestCase

from light.landmarks import ALL_LANDMARK_FINDER_CLASSES
from light.trig import ALL_TRIG_FINDER_CLASSES


class TestSymbols(TestCase):
    """
    Tests to ensure that all landmarks and trig points have different symbols.
    """
    def testRightNumberOfSymbols(self):
        """
        Must find the right number of symbols.
        """
        trigSymbols = [symbol.SYMBOL for symbol in ALL_TRIG_FINDER_CLASSES]
        lmSymbols = [symbol.SYMBOL for symbol in ALL_LANDMARK_FINDER_CLASSES]
        symbols = trigSymbols + lmSymbols
        self.assertEqual(10, len(symbols))

    def testAllLandmarkSymbolsAreDifferent(self):
        """
        All landmark symbols must be different.
        """
        symbols = [symbol.SYMBOL for symbol in ALL_LANDMARK_FINDER_CLASSES]
        self.assertEqual(len(symbols), len(set(symbols)))

    def testAllSymbolsAreDifferent(self):
        """
        All trig finder symbols must be different.
        """
        symbols = [symbol.SYMBOL for symbol in ALL_TRIG_FINDER_CLASSES]
        self.assertEqual(len(symbols), len(set(symbols)))
