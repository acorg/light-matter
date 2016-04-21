from unittest import TestCase

from light.performance.polymerase import BIT_SCORES, Z_SCORES

POLYMERASE_COUNT = 22

POLYMERASE_NAMES = {
    '2J7W-A', '4K6M-A', '1S49-A', '1NB6-A', '3OLB-A', '1XR7-A', '1XR6-A',
    '3CDW-A', '2E9Z-A', '3BSO-A', '3UQS-A', '1KHV-B', '2CKW-A', '1HI0-P',
    '3AVX-A', '2PUS-A', '2YI9-A', '2R7W-A', '1N35-A', '3V81-C', '1MU2-A',
    '1RW3-A'
}


class TestZScores(TestCase):
    """
    Tests for the light.performance.polymerase Z_SCORES variable.
    """
    def testVariables(self):
        """
        Make sure the variables defined above are consistent.
        """
        self.assertEqual(POLYMERASE_COUNT, len(POLYMERASE_NAMES))

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
                ('3BSO-A', '3OLB-A', 32.0),
                ('2J7W-A', '4K6M-A', 42.9),
                ('1RW3-A', '1MU2-A', 20.7),
                ('2YI9-A', '1HI0-P', 10.7)):
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
        self.assertEqual(POLYMERASE_NAMES, set(BIT_SCORES.keys()))

    def testValues(self):
        """
        BIT_SCORES values must have the correct keys. Each key, k (a query id),
        must have a set of keys (subject ids) that is the full set of sequence
        ids except for k.
        """
        for queryId in POLYMERASE_NAMES:
            expected = POLYMERASE_NAMES - {queryId}
            self.assertEqual(expected, set(BIT_SCORES[queryId]))

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
