from unittest import TestCase

from performance.bin.evaluate_helices import evaluateMatch


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
        self.assertFalse(evaluateMatch('GSSSGGGGG', 1, 7))

    def testCase2(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence doesn't extend to the left or right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('SHHHHHHG', 1, 7))

    def testCase3(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to the left.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('HHHHHHHG', 1, 7))

    def testCase4(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to the right.
        --> true positive.
        """
        self.assertTrue(evaluateMatch('GHHHHHHH', 1, 7))

    def testHelixExtendsBothSites(self):
        """
        The alpha helix matches part of a sequence that's an alpha helix. The
        alpha helix in the sequence extends to both sides.
        --> false positive.
        """
        self.assertFalse(evaluateMatch('HHHHHHHH', 1, 7))