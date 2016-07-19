from unittest import TestCase

from light.performance.evaluate import evaluateMatch, evaluateMatchNoPrefix


class TestEvaluateMatch(TestCase):
    """
    Tests for the evaluateMatch function. For a description of the cases tested
    see the docstring of evaluateMatch.
    """
    def testCase1(self):
        """
        The alpha helix matches part of a sequence that's not an alpha helix.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('GSSSGGGGG', 1, 7, 'H'))

    def testCase1NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that's not that structure.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('GSSSGGGGG', 1, 7, 'E'))

    def testCase2(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence doesn't extend to the left or right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('SHHHHHHG', 1, 7, 'H'))

    def testCase2NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that isn't that structure. The
        structure in the sequence doesn't extend to the left or right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('SEEEEEEG', 1, 7, 'E'))

    def testCase2ExactMatch(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence doesn't extend to the left or right.
        --> true positive.
        The helix starts at 0 and matches exactly.
        """
        self.assertTrue(evaluateMatch('HHH', 0, 3, 'H'))

    def testCase3(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to the left.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('HHHHHHHG', 1, 7, 'H'))

    def testCase3NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that's that structure. The
        structure in the sequence extends to the left.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('EEEEEEEG', 1, 7, 'E'))

    def testCase4(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to the right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('GHHHHHHH', 1, 7, 'H'))

    def testCase4NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that's that structure. The
        structure in the sequence extends to the right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('GEEEEEEE', 1, 7, 'E'))

    def testCase4HelixStartsAt0(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence doesn't extends to the right.
        --> true positive.
        The helix starts at 0 and extends to the right.
        """
        self.assertTrue(evaluateMatch('HHHHHHHH', 0, 6, 'H'))

    def testHelixExtendsBothSides(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to both sides.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('HHHHHHHH', 1, 7, 'H'))

    def testAlphaHelixCombined(self):
        """
        The alpha helix combined matches part of a sequence that's an alpha
        helix.
        """
        self.assertTrue(evaluateMatch('HHGHGHIH', 1, 7, 'K'))


class TestEvaluateMatchNoPrefix(TestCase):
    """
    Tests for the evaluateMatchNoPrefix function. For a description of the
    cases tested see the docstring of evaluateMatchNoPrefix.
    """
    def testCase1(self):
        """
        The alpha helix matches part of a sequence that's not an alpha helix.
        --> false positive.
        """
        self.assertFalse(evaluateMatchNoPrefix('GSSSGGGGG', 1, 7, 'H'))

    def testCase1NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that's not that structure.
        --> false positive.
        """
        self.assertFalse(evaluateMatchNoPrefix('ISSSIIIII', 1, 7, 'E'))

    def testCase2(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix.
        --> true positive.
        """
        self.assertTrue(evaluateMatchNoPrefix('SHHHHHHG', 1, 7, 'H'))

    def testCase2NotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that isn't that structure.
        --> true positive.
        """
        self.assertTrue(evaluateMatchNoPrefix('SEEEEEEG', 1, 7, 'E'))

    def testCase2ExactMatch(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix.
        --> true positive.
        The helix starts at 0 and matches exactly.
        """
        self.assertTrue(evaluateMatchNoPrefix('HHH', 0, 3, 'H'))

    def testCase2ExactMatchNotDefault(self):
        """
        Test that the structure evaluation works when a non-default
        structureType is used.
        The structure matches part of a sequence that isn't that structure.
        --> true positive.
        The structure starts at 0 and matches exactly.
        """
        self.assertTrue(evaluateMatchNoPrefix('EEE', 0, 3, 'E'))

    def testAlphaHelixCombinedMustBeFound(self):
        """
        The alpha helix combined matches part of a sequence that's an alpha
        helix.
        """
        self.assertTrue(evaluateMatchNoPrefix('GHHHGGIGGG', 1, 7, 'K'))
