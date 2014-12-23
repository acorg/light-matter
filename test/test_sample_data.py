from json import loads
from unittest import TestCase

from dark.hsp import HSP

from light.alignments import jsonDictToAlignments

from sample_data import (
    DB, PARAMS, COWPOX, MONKEYPOX, MUMMYPOX, SQUIRRELPOX55, SQUIRRELPOX1296,
    RECORD0, RECORD1, RECORD2, RECORD3, RECORD4)


class TestParams(TestCase):
    """
    Test the database parameters.
    """
    def testPARAMS(self):
        """
        The saved (JSON format) database parameters must be as expected.
        """
        self.maxDiff = None
        self.assertEqual(
            {
                'subjectCount': 5,
                'totalCoveredResidues': 61,
                'checksum': ('e00eba538693d82ef6d275cc935c5b55'
                             'aede072c040fe312a47bbc268ba10f56'),
                'limitPerLandmark': None,
                'trigPointFinderClasses': [
                    'AminoAcids',
                    'Troughs',
                ],
                'totalResidues': 68,
                'landmarkFinderClasses': [
                    'AlphaHelix',
                    'AlphaHelix_3_10',
                    'AlphaHelix_pi',
                    'BetaStrand',
                ],
                'maxDistance': None,
            },
            loads(PARAMS))


class TestRecords(TestCase):
    """
    Test the result of running find on the reads in the sample data.
    """
    def testRECORD0(self):
        """
        The results in RECORD0 must match two subjects: SQUIRRELPOX1296 (with
        score 2) and SQUIRRELPOX55 (with score 2).
        """
        alignments = jsonDictToAlignments(loads(RECORD0), DB)
        self.assertEqual(2, len(alignments))

        self.assertEqual(SQUIRRELPOX1296.id, alignments[0].subjectTitle)
        self.assertEqual(HSP(5), alignments[0].hsps[0])

        self.assertEqual(SQUIRRELPOX55.id, alignments[1].subjectTitle)
        self.assertEqual(HSP(4), alignments[1].hsps[0])

    def testRECORD1(self):
        """
        The results in RECORD1 must match two subjects: MONKEYPOX (with score
        1) and MUMMYPOX (with score 1).
        """
        alignments = jsonDictToAlignments(loads(RECORD1), DB)
        self.assertEqual(2, len(alignments))

        self.assertEqual(MONKEYPOX.id, alignments[0].subjectTitle)
        self.assertEqual(HSP(3), alignments[0].hsps[0])

        self.assertEqual(MUMMYPOX.id, alignments[1].subjectTitle)
        self.assertEqual(HSP(1), alignments[1].hsps[0])

    def testRECORD2(self):
        """
        The results in RECORD2 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD2), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(HSP(36), alignments[0].hsps[0])

    def testRECORD3(self):
        """
        The results in RECORD3 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD3), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(HSP(36), alignments[0].hsps[0])

    def testRECORD4(self):
        """
        The results in RECORD4 must match just COWPOX.
        """
        alignments = jsonDictToAlignments(loads(RECORD4), DB)
        self.assertEqual(1, len(alignments))

        self.assertEqual(COWPOX.id, alignments[0].subjectTitle)
        self.assertEqual(HSP(36), alignments[0].hsps[0])
