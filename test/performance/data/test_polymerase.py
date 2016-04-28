from unittest import TestCase

from light.performance.data.polymerase import (
    QUERIES, SUBJECTS, BIT_SCORES, Z_SCORES)

QUERY_COUNT = SUBJECT_COUNT = 22

QUERY_NAMES = {
    'pdb_2j7w_a', 'pdb_4k6m_a', 'pdb_1s49_a', 'pdb_1nb6_a', 'pdb_3olb_a',
    'pdb_1xr7_a', 'pdb_1xr6_a', 'pdb_3cdw_a', 'pdb_2e9z_a', 'pdb_3bso_a',
    'pdb_3uqs_a', 'pdb_1khv_b', 'pdb_2ckw_a', 'pdb_1hi0_p', 'pdb_3avx_a',
    'pdb_2pus_a', 'pdb_2yi9_a', 'pdb_2r7w_a', 'pdb_1n35_a', 'pdb_3v81_c',
    'pdb_1mu2_a', 'pdb_1rw3_a'
}


class TestQueries(TestCase):
    """
    Tests for the light.performance.pdb_2hla_a QUERIES variable.
    """
    def testLen(self):
        """
        Make sure the query counts are as expected.
        """
        self.assertEqual(QUERY_COUNT, len(QUERIES))

    def testNoDuplicates(self):
        """
        There shouldn't be any duplicated queries. If the query has known
        structure, both the sequence and structure are part of what's
        matched in detecting duplicates. This allows identical sequences with
        different structures to be part of a query set.
        """
        seen = set()
        for query in QUERIES:
            s = query.sequence + '\0' + getattr(query, 'structure', '')
            self.assertNotIn(s, seen)
            seen.add(s)


class TestSubjects(TestCase):
    """
    Tests for the light.performance.pdb_2hla_a SUBJECTS variable.

    In theory we don't need to do these tests since for this data set the
    queries and subjects are identical. But let's be defensive.
    """
    def testLen(self):
        """
        Make sure the subject counts are as expected.
        """
        self.assertEqual(SUBJECT_COUNT, len(SUBJECTS))

    def testNoDuplicates(self):
        """
        There shouldn't be any duplicated subjects. If the subject has known
        structure, both the sequence and structure are part of what's
        matched in detecting duplicates. This allows identical sequences with
        different structures to be part of a subject set.
        """
        seen = set()
        for subject in SUBJECTS:
            s = subject.sequence + '\0' + getattr(subject, 'structure', '')
            self.assertNotIn(s, seen)
            seen.add(s)


class TestZScores(TestCase):
    """
    Tests for the light.performance.data.polymerase Z_SCORES variable.
    """
    def testVariables(self):
        """
        Make sure the variables defined above are consistent.
        """
        self.assertEqual(QUERY_COUNT, len(QUERY_NAMES))

    def testLen(self):
        """
        Z_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(QUERY_COUNT, len(Z_SCORES))

        for queryId in Z_SCORES:
            self.assertEqual(QUERY_COUNT - 1, len(Z_SCORES[queryId]))

    def testLogic(self):
        """
        The Z_SCORES dict must contain the expected values.

        Here we make reasonably sure that the rearrangement of _ZSCORE_TABLE in
        light/performance/polymerase.py is correct. This is just a selection
        of elements chosen from the overall table to make sure things look
        right.
        """
        for queryId, subjectId, expected in (
                ('pdb_3bso_a', 'pdb_3olb_a', 32.0),
                ('pdb_2j7w_a', 'pdb_4k6m_a', 42.9),
                ('pdb_1rw3_a', 'pdb_1mu2_a', 20.7),
                ('pdb_2yi9_a', 'pdb_1hi0_p', 10.7)):
            self.assertAlmostEqual(expected, Z_SCORES[queryId][subjectId])

    def testSymmetric(self):
        """
        Z_SCORES must be symmetric.
        """
        for queryId in QUERY_NAMES:
            for subjectId in QUERY_NAMES:
                if queryId != subjectId:
                    self.assertEqual(Z_SCORES[queryId][subjectId],
                                     Z_SCORES[subjectId][queryId])

    def testSelfScoresAbsent(self):
        """
        Z_SCORES must not exist for sequences matched against themselves.
        """
        for queryId in QUERY_NAMES:
            self.assertNotIn(queryId, Z_SCORES[queryId])

    def testNonNegative(self):
        """
        Z_SCORES cannot be less than zero.
        """
        for queryId in QUERY_NAMES:
            for subjectId in QUERY_NAMES:
                if queryId != subjectId:
                    self.assertTrue(Z_SCORES[queryId][subjectId] >= 0.0)


class TestBitScores(TestCase):
    """
    Tests for the light.performance.data.polymerase BIT_SCORES variable.
    """
    def testLen(self):
        """
        BIT_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(QUERY_COUNT, len(BIT_SCORES))

        for queryId in BIT_SCORES:
            self.assertEqual(QUERY_COUNT, len(BIT_SCORES[queryId]))

    def testKeys(self):
        """
        BIT_SCORES must have the correct keys.
        """
        self.assertEqual(QUERY_NAMES, set(BIT_SCORES))

    def testValues(self):
        """
        BIT_SCORES values must have the correct keys. Each key, k (a query id),
        must have a set of keys (subject ids) that is the full set of sequence
        ids.
        """
        for queryId in QUERY_NAMES:
            self.assertEqual(QUERY_NAMES, set(BIT_SCORES[queryId]))

    def testSymmetric(self):
        """
        BIT_SCORES must be symmetric.
        """
        for queryId in QUERY_NAMES:
            for subjectId in QUERY_NAMES:
                if queryId != subjectId:
                    self.assertEqual(BIT_SCORES[queryId][subjectId],
                                     BIT_SCORES[subjectId][queryId])

    def testNonNegative(self):
        """
        BIT_SCORES cannot be less than zero.
        """
        for queryId in QUERY_NAMES:
            for subjectId in QUERY_NAMES:
                if queryId != subjectId:
                    self.assertTrue(BIT_SCORES[queryId][subjectId] >= 0.0)
