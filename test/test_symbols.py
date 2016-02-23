from unittest import TestCase

from light.landmarks import ALL_LANDMARK_CLASSES_INCLUDING_DEV
from light.trig import ALL_TRIG_CLASSES_INCLUDING_DEV


class TestSymbols(TestCase):
    """
    Tests to ensure that all landmarks and trig points have different symbols.
    """
    def testRightNumberOfSymbols(self):
        """
        Must find the right number of symbols.
        """
        symbols = [cls.SYMBOL for cls in
                   ALL_LANDMARK_CLASSES_INCLUDING_DEV +
                   ALL_TRIG_CLASSES_INCLUDING_DEV]
        self.assertEqual(22, len(symbols))

    def testAllSymbolsAreDistinct(self):
        """
        All landmark and trig finder symbols must be distinct.
        """
        symbols = [cls.SYMBOL for cls in
                   ALL_LANDMARK_CLASSES_INCLUDING_DEV +
                   ALL_TRIG_CLASSES_INCLUDING_DEV]
        self.assertEqual(len(symbols), len(set(symbols)))
