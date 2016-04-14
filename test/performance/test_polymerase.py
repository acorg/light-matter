from unittest import TestCase

from light.performance.polymerase import BIT_SCORES, Z_SCORES

POLYMERASE_COUNT = 22

# Note that the names must be in this order. If they're not, the code in
# light/performance/polymerase.py is incorrect.
POLYMERASE_NAMES = [
    '2J7W', '4K6M', '1S49', '1NB6', '3OLB', '1XR7', '1XR6', '3CDW', '2E9Z',
    '3BSO', '3UQS', '1KHV', '2CKW', '1HI0', '3AVX', '2PUS', '2YI9', '2R7W',
    '1N35', '3V81', '1MU2', '1RW3'
]


class TestZScores(TestCase):
    """
    Tests for the light.performance.polymerase Z_SCORES variable.
    """
    def testLen(self):
        """
        Z_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(POLYMERASE_COUNT, len(Z_SCORES))

        for queryId in Z_SCORES:
            self.assertEqual(POLYMERASE_COUNT - 1, len(Z_SCORES[queryId]))

    def testLogic(self):
        """
        The Z_SCORES dict must contain the expected values.

        Here we make reasonably sure that the rearrangement of _ZSCORE_TABLE in
        light/performance/polymerase.py is correct. This is just a selection
        of elements chosen from the overall table to make sure things look
        right.
        """
        for queryId, subjectId, expected in (
                ('3BSO', '3OLB', 32.0),
                ('2J7W', '4K6M', 42.9),
                ('1RW3', '1MU2', 20.7),
                ('2YI9', '1HI0', 10.7)):
            self.assertAlmostEqual(expected, Z_SCORES[queryId][subjectId])

    def testSymmetric(self):
        """
        Z_SCORES must be symmetric.
        """
        for queryId in POLYMERASE_NAMES:
            for subjectId in POLYMERASE_NAMES:
                if queryId != subjectId:
                    self.assertEqual(Z_SCORES[queryId][subjectId],
                                     Z_SCORES[subjectId][queryId])

    def testSelfScoresAbsent(self):
        """
        Z_SCORES must not exist for sequences matched against themselves.
        """
        for queryId in POLYMERASE_NAMES:
            self.assertNotIn(queryId, Z_SCORES[queryId])

    def testNonNegative(self):
        """
        Z_SCORES cannot be less than zero.
        """
        for queryId in POLYMERASE_NAMES:
            for subjectId in POLYMERASE_NAMES:
                if queryId != subjectId:
                    self.assertTrue(Z_SCORES[queryId][subjectId] >= 0.0)


class TestBitScores(TestCase):
    """
    Tests for the light.performance.polymerase BIT_SCORES variable.
    """
    def testLen(self):
        """
        BIT_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(POLYMERASE_COUNT, len(BIT_SCORES))

        for queryId in BIT_SCORES:
            self.assertEqual(POLYMERASE_COUNT - 1, len(BIT_SCORES[queryId]))

    def testKeys(self):
        """
        BIT_SCORES must have the correct keys.
        """
        self.assertEqual(set(POLYMERASE_NAMES), set(BIT_SCORES.keys()))

    def testValues(self):
        """
        BIT_SCORES values must have the correct keys.
        """
        for queryId in POLYMERASE_NAMES:
            expected = set(POLYMERASE_NAMES)
            expected.remove(queryId)
            self.assertEqual(set(expected), set(BIT_SCORES[queryId]))

    def testSymmetric(self):
        """
        BIT_SCORES must be symmetric.
        """
        for queryId in POLYMERASE_NAMES:
            for subjectId in POLYMERASE_NAMES:
                if queryId != subjectId:
                    self.assertEqual(BIT_SCORES[queryId][subjectId],
                                     BIT_SCORES[subjectId][queryId])

    def testSelfScoresAbsent(self):
        """
        BIT_SCORES must not exist for sequences matched against themselves.
        """
        for queryId in POLYMERASE_NAMES:
            self.assertNotIn(queryId, BIT_SCORES[queryId])

    def testNonNegative(self):
        """
        BIT_SCORES cannot be less than zero.
        """
        for queryId in POLYMERASE_NAMES:
            for subjectId in POLYMERASE_NAMES:
                if queryId != subjectId:
                    self.assertTrue(BIT_SCORES[queryId][subjectId] >= 0.0)
