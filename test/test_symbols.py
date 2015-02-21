from unittest import TestCase

from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES


class TestSymbols(TestCase):
    """
    Tests to ensure that all landmarks and trig points have different symbols.
    """
    def testRightNumberOfSymbols(self):
        """
        Must find the right number of symbols.
        """
        symbols = ([cls.SYMBOL for cls in ALL_LANDMARK_CLASSES] +
                   [cls.SYMBOL for cls in ALL_TRIG_CLASSES])
        self.assertEqual(15, len(symbols))

    def testAllSymbolsAreDifferent(self):
        """
        All landmark and trig finder symbols must be different.
        """
        symbols = ([cls.SYMBOL for cls in ALL_LANDMARK_CLASSES] +
                   [cls.SYMBOL for cls in ALL_TRIG_CLASSES])
        self.assertEqual(len(symbols), len(set(symbols)))
