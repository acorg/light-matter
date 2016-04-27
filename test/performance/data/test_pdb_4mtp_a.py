from unittest import TestCase, skip

from light.performance.data.pdb_4mtp_a import (
    QUERIES, SUBJECTS, BIT_SCORES, Z_SCORES)

SUBJECT_COUNT = 1
QUERY_COUNT = 217


class TestQueries(TestCase):
    """
    Tests for the light.performance.pdb_2hla_a QUERIES variable.
    """
    def testLen(self):
        """
        Make sure the query counts are as expected.
        """
        self.assertEqual(QUERY_COUNT, len(QUERIES))


class TestSubjects(TestCase):
    """
    Tests for the light.performance.pdb_2hla_a SUBJECTS variable.
    """
    def testLen(self):
        """
        Make sure the subject counts are as expected.
        """
        self.assertEqual(SUBJECT_COUNT, len(SUBJECTS))


class TestBitScores(TestCase):
    """
    Tests for the light.performance.data.pdb_2hla_a BIT_SCORES variable.
    """
    def testLen(self):
        """
        BIT_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(QUERY_COUNT, len(BIT_SCORES))

        for queryId in BIT_SCORES:
            self.assertEqual(SUBJECT_COUNT, len(BIT_SCORES[queryId]))

    def testKeys(self):
        """
        BIT_SCORES must have the correct keys.
        """
        self.assertEqual(set(query.id for query in QUERIES), set(BIT_SCORES))

    def testValues(self):
        """
        BIT_SCORES values must have the correct keys. Each key, k (a query id),
        must have a set of keys (subject ids) that is the full set of sequence
        ids.
        """
        subjectIds = set(subject.id for subject in SUBJECTS)
        for queryId in BIT_SCORES:
            self.assertEqual(subjectIds, set(BIT_SCORES[queryId]))

    def testNonNegative(self):
        """
        BIT_SCORES cannot be less than zero.
        """
        subjectIds = set(subject.id for subject in SUBJECTS)
        for queryId in BIT_SCORES:
            for subjectId in subjectIds:
                if queryId != subjectId:
                    self.assertTrue(BIT_SCORES[queryId][subjectId] >= 0.0)


class TestZScores(TestCase):
    """
    Tests for the light.performance.pdb_2hla_a Z_SCORES variable.
    """
    @skip('Pending investigation into missing FASTA seqs')
    def testLen(self):
        """
        Z_SCORES and each of its keys must have the expected number of
        values.
        """
        self.assertEqual(QUERY_COUNT, len(Z_SCORES))

        for queryId in Z_SCORES:
            self.assertEqual(SUBJECT_COUNT, len(Z_SCORES[queryId]))

    @skip('Pending investigation into missing FASTA seqs')
    def testKeys(self):
        """
        Z_SCORES must have the correct keys.
        """
        self.assertEqual(set(query.id for query in QUERIES), set(Z_SCORES))

    def testValues(self):
        """
        Z_SCORES values must have the correct keys. Each key, k (a query id),
        must have a set of keys (subject ids) that is the full set of sequence
        ids.
        """
        subjectIds = set(subject.id for subject in SUBJECTS)
        for queryId in Z_SCORES:
            self.assertEqual(subjectIds, set(Z_SCORES[queryId]))

    def testNonNegative(self):
        """
        Z_SCORES cannot be less than zero.
        """
        for queryId in Z_SCORES:
            for subjectId in Z_SCORES[queryId]:
                self.assertTrue(Z_SCORES[queryId][subjectId] >= 0.0)
