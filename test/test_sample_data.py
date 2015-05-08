from json import loads
from unittest import TestCase

from light.alignments import jsonDictToAlignments

from .sample_data import (
    DB, COWPOX, MONKEYPOX, MUMMYPOX, CATPOX, SQUIRRELPOX,
    RECORD0, RECORD1, RECORD2, RECORD3, RECORD4,
    READ0_SQUIRRELPOX_SCORE, READ0_CATPOX_SCORE,
    READ1_MONKEYPOX_SCORE, READ1_MONKEYPOX_HSP2_SCORE, READ1_MUMMYPOX_SCORE,
    READ2_COWPOX_SCORE,
    READ3_COWPOX_SCORE,
    READ4_COWPOX_SCORE)


class TestRecords(TestCase):
    """
    Test the result of running find on the reads in the sample data.
    """
    def testRECORD0(self):
        """
        The results in RECORD0 must match two subjects: CATPOX
        and SQUIRRELPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD0), DB)
        self.assertEqual(2, len(alignments))

        self.assertEqual(CATPOX.id, alignments[0].subjectTitle)
        self.assertEqual(READ0_CATPOX_SCORE, alignments[0].hsps[0].score.score)

        self.assertEqual(SQUIRRELPOX.id, alignments[1].subjectTitle)
        self.assertEqual(READ0_SQUIRRELPOX_SCORE,
                         alignments[1].hsps[0].score.score)

    def testRECORD1(self):
        """
        The results in RECORD1 must match two subjects: MONKEYPOX and MUMMYPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD1), DB)
        self.assertEqual(2, len(alignments))

        self.assertEqual(MONKEYPOX.id, alignments[0].subjectTitle)
        self.assertEqual(READ1_MONKEYPOX_SCORE,
                         alignments[0].hsps[0].score.score)

        self.assertEqual(MUMMYPOX.id, alignments[1].subjectTitle)
        self.assertEqual(READ1_MUMMYPOX_SCORE,
                         alignments[1].hsps[0].score.score)

    def testRECORD2(self):
        """
        The results in RECORD2 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD2), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(READ2_COWPOX_SCORE, alignments[0].hsps[0].score.score)

    def testRECORD3(self):
        """
        The results in RECORD3 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD3), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(READ3_COWPOX_SCORE, alignments[0].hsps[0].score.score)

    def testRECORD4(self):
        """
        The results in RECORD4 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD4), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(READ4_COWPOX_SCORE, alignments[0].hsps[0].score.score)


class TestSubjects(TestCase):
    """
    Test the subject sequences in the sample database. These tests are not
    needed for correctness. They are included to show the specific values
    of sequence lengths and scores and so that we will have failing tests
    if any of those things change unexpectedly.

    Having the explicit numerical values also helps to see why tests elsewhere
    have the results they do.
    """
    def testCATPOX(self):
        """
        The CATPOX subject must have the expected length and match scores.
        """
        self.assertEqual(34, len(CATPOX.sequence))
        self.assertEqual(0.5, READ0_CATPOX_SCORE)

    def testCOWPOX(self):
        """
        The COWPOX subject must have the expected length and match scores.
        """
        self.assertEqual(14, len(COWPOX.sequence))
        self.assertEqual(0.75, READ2_COWPOX_SCORE)
        self.assertEqual(0.75, READ3_COWPOX_SCORE)
        self.assertEqual(0.75, READ4_COWPOX_SCORE)

    def testMONKEYPOX(self):
        """
        The MONKEYPOX subject must have the expected length and match scores.
        """
        self.assertEqual(47, len(MONKEYPOX.sequence))
        self.assertEqual(0.6, READ1_MONKEYPOX_SCORE)
        self.assertEqual(0.3, READ1_MONKEYPOX_HSP2_SCORE)

    def testMUMMYPOX(self):
        """
        The MUMMYPOX subject must have the expected length and match scores.
        """
        self.assertEqual(48, len(MUMMYPOX.sequence))
        self.assertEqual(0.4, READ1_MUMMYPOX_SCORE)

    def testSQUIRRELPOX(self):
        """
        The SQUIRRELPOX subject must have the expected length and match scores.
        """
        self.assertEqual(41, len(SQUIRRELPOX.sequence))
        self.assertEqual(1.0, READ0_SQUIRRELPOX_SCORE)
