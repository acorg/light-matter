import warnings
from unittest import TestCase
import numpy as np

from light.parameters import FindParameters
from light.graphics import PlotHashesInSubjectAndRead, SequenceFeatureAnalysis

from dark.reads import AARead


class TestPlotHashesInSubjectAndRead(TestCase):
    """
    Tests for the light.graphics.PlotHashesInSubjectAndRead class.
    """
    def testNoHashes(self):
        """
        If there are no hashes in subject and query, matchingHashes,
        queryHashes and subjectHashes must all be empty.
        """
        subject = AARead('subject', 'RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR')
        query = AARead('query', 'AAAAAAAAAAAAAAAAAAAAAAA')
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_3_10'],
            trigPointNames=['Peaks'])
        self.assertEqual(0, len(list(hashes.matchingHashes)))
        self.assertEqual(0, len(list(hashes.queryHashes)))
        self.assertEqual(0, len(list(hashes.subjectHashes)))

    def testNoMatchingHashes(self):
        """
        If there are no matching hashes in subject and query, matchingHashes
        must be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAAASARRRRRRRRRRRRRR')
        query = AARead('query', 'FRRFRRFAAASAAAAAAAAAAAAA')
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_3_10'],
            trigPointNames=['Peaks'])
        self.assertEqual(0, len(list(hashes.matchingHashes)))
        self.assertEqual(1, len(list(hashes.queryHashes)))
        self.assertEqual(1, len(list(hashes.subjectHashes)))

    def testNoQueryOrSubjectHashes(self):
        """
        If all hashes in the query also occur in the subject, queryHashes must
        be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAASARRRRFRRFRRFAAASA')
        query = AARead('query', 'ASARRRRFRRFRRFAAASA')
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_3_10'],
            trigPointNames=['Peaks'])
        self.assertEqual(2, len(hashes.matchingHashes))
        self.assertEqual(0, len(hashes.queryHashes))
        self.assertEqual(3, len(hashes.subjectHashes))

    def testNoSubjectHashes(self):
        """
        If all hashes in the subject also occur in the query, subjectHashes
        must be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAASA')
        query = AARead('query', 'FRRFRRFRRFAAAAAAAAAAAASARRRRFRRFRRFAAASA')
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_3_10'],
            trigPointNames=['Peaks'])
        self.assertEqual(1, len(hashes.matchingHashes))
        self.assertEqual(4, len(hashes.queryHashes))
        self.assertEqual(0, len(hashes.subjectHashes))

    def testShowSignificant(self):
        """
        The showSignificant option must work correctly.
        """
        seq = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')

        # Showing the significant:
        hashes = PlotHashesInSubjectAndRead(
            seq, seq, landmarkNames=['AlphaHelix', 'BetaStrand'],
            trigPointNames=['Peaks'])
        self.assertEqual(13, len(hashes.matchingHashes))
        self.assertEqual(0, len(hashes.queryHashes))
        self.assertEqual(0, len(hashes.subjectHashes))

        # Same input, but hiding the significant:
        hashes = PlotHashesInSubjectAndRead(
            seq, seq, landmarkNames=['AlphaHelix', 'BetaStrand'],
            trigPointNames=['Peaks'], showSignificant=False)
        self.assertEqual(1, len(hashes.matchingHashes))
        self.assertEqual(0, len(hashes.queryHashes))
        self.assertEqual(12, len(hashes.subjectHashes))

    def testShowInsignificant(self):
        """
        The showInsignificant option must work correctly.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASAVVVVVVASAVVVASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')

        # Showing the insignificant:
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'BetaStrand'],
            trigPointNames=['Peaks'])
        self.assertEqual(2, len(hashes.matchingHashes))
        self.assertEqual(11, len(hashes.queryHashes))
        self.assertEqual(7, len(hashes.subjectHashes))

        # Same input, but hiding the significant:
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'BetaStrand'],
            trigPointNames=['Peaks'], showInsignificant=False)
        self.assertEqual(0, len(hashes.matchingHashes))
        self.assertEqual(11, len(hashes.queryHashes))
        self.assertEqual(9, len(hashes.subjectHashes))

    def testShowBestBinOnly(self):
        """
        The showBestBinOnly option must work correctly, by setting the 'bins'
        attribute to contain just one bin.
        """

        A = 'FRRRFRRRFXXXXXX'
        C = 'FRRRRFRRRRFXXXXXX'

        subject = AARead('subject', 5 * A + C + 2 * A)
        query = AARead('query', 5 * A)

        findParams = FindParameters(significanceFraction=0.01)

        # There are 11 significant bins.
        hashes = PlotHashesInSubjectAndRead(
            query, subject, findParams,
            landmarkNames=['AlphaHelix', 'AlphaHelix_pi'], trigPointNames=[],
            distanceBase=1.025, limitPerLandmark=50, minDistance=1,
            maxDistance=100, showInsignificant=False)
        self.assertEqual(11, len(hashes.bins))

        bestScore = hashes.result.analysis['bestScore']

        # Same input, but restricting ourselves to only the single most
        # significant bin:
        hashes = PlotHashesInSubjectAndRead(
            query, subject, findParams,
            landmarkNames=['AlphaHelix', 'AlphaHelix_pi'], trigPointNames=[],
            distanceBase=1.025, limitPerLandmark=50, minDistance=1,
            maxDistance=100, showInsignificant=False, showBestBinOnly=True)
        self.assertEqual(1, len(hashes.bins))

        # Check that the best bin when we use onlyShowBestBin has the same
        # score as the bin we get when we don't use onlyShowBestBin.
        self.assertEqual(bestScore, hashes.result.analysis['bestScore'])

    def testShowBestBinOnlyIssuesWarning(self):
        """
        The showBestBinOnly option must issue a warning when more than one bin
        has the best score.
        """

        A = 'FRRRFRRRFXXXXXX'
        C = 'FRRRRFRRRRFXXXXXX'

        subject = AARead('subject', 5 * A + C + 5 * A)
        query = AARead('query', 5 * A)

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore')

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            PlotHashesInSubjectAndRead(
                query, subject, findParams,
                landmarkNames=['AlphaHelix', 'AlphaHelix_pi'],
                trigPointNames=[], distanceBase=1.025, limitPerLandmark=50,
                minDistance=1, maxDistance=100, showInsignificant=False,
                showBestBinOnly=True)

            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Multiple bins share the best score (1.000000). '
                     'Displaying just one of them.')
            self.assertIn(error, str(w[0].message))


class TestSequenceFeatureAnalysis(TestCase):
    """
    Tests for the light.graphics.SequenceFeatureAnalysis class.
    """
    def testOneFeatureUniqueOffsets(self):
        """
        If the sequence contains a single feature, the unique offsets method
        must return the correct value.
        """
        read = AARead('id', 'FRRRFRRRF')
        s = SequenceFeatureAnalysis(read, landmarkNames=['AlphaHelix'],
                                    trigPointNames=['Peaks'])
        self.assertEqual(
            {
                'AlphaHelix': {0, 1, 2, 3, 4, 5, 6, 7, 8},
                'Peaks': set(),
            },
            s.uniqueOffsets())

    def testTwoNonOverlappingFeatureUniqueOffsets(self):
        """
        If the sequence contains two non-overlapping features, the unique
        offsets method must return the correct value.
        """
        read = AARead('id', 'FRRRFRRRFW')
        s = SequenceFeatureAnalysis(read, landmarkNames=['AlphaHelix'],
                                    trigPointNames=['AminoAcids'])
        self.assertEqual(
            {
                'AlphaHelix': {0, 1, 2, 3, 4, 5, 6, 7, 8},
                'AminoAcids': {9},
            },
            s.uniqueOffsets())

    def testTwoSpreadOutNonOverlappingFeatureUniqueOffsets(self):
        """
        If the sequence contains two non-overlapping features with other AAs
        before and after them, the unique offsets method must return the
        correct value.
        """
        read = AARead('id', 'XXFRRRFRRRFXXWX')
        s = SequenceFeatureAnalysis(read, landmarkNames=['AlphaHelix'],
                                    trigPointNames=['AminoAcids'])
        self.assertEqual(
            {
                'AlphaHelix': {2, 3, 4, 5, 6, 7, 8, 9, 10},
                'AminoAcids': {13},
            },
            s.uniqueOffsets())

    def testTwoOverlappingFeatureUniqueOffsets(self):
        """
        If the sequence contains two features that overlap, the unique offsets
        method must return the correct value, which will exclude the offset(s)
        in common.
        """
        read = AARead('id', 'AGYGSTWT')
        s = SequenceFeatureAnalysis(read,
                                    landmarkNames=['AlphaHelix', 'Prosite'],
                                    trigPointNames=['AminoAcids'])
        # There is an AminoAcid trig point at offset 6, but that's also in
        # the Prosite feature, so offset 6 is not unique to either finder.
        self.assertEqual(
            {
                'AlphaHelix': set(),
                'AminoAcids': set(),
                'Prosite': {0, 1, 2, 3, 4, 5, 7},
            },
            s.uniqueOffsets())

    def testNames(self):
        """
        The __init__ method of SequenceFeatureAnalysis must set the correct
        landmark and trig point names based on what's in the sequence and what
        we were looking for.
        """
        read = AARead('id', 'FRRRFRRRFWNPNW')
        s = SequenceFeatureAnalysis(read,
                                    landmarkNames=['AlphaHelix', 'BetaTurn'],
                                    trigPointNames=['AminoAcids'])
        self.assertEqual(set(['AlphaHelix', 'BetaTurn']), s.landmarkNames)
        self.assertEqual(set(['AminoAcids']), s.trigPointNames)

    def testMissingNames(self):
        """
        The __init__ method of SequenceFeatureAnalysis must set the correct
        missing landmark and trig point names based on what's in the sequence
        and what we were looking for.
        """
        read = AARead('id', 'FRRRFRRRF')
        s = SequenceFeatureAnalysis(read,
                                    landmarkNames=['AlphaHelix', 'BetaTurn'],
                                    trigPointNames=['AminoAcids'])
        self.assertEqual(set(['BetaTurn']), s.landmarksNotFound)
        self.assertEqual(set(['AminoAcids']), s.trigPointsNotFound)

    def testOffsets(self):
        """
        The __init__ method of SequenceFeatureAnalysis must set the correct
        landmark and trig point offsets.
        """
        read = AARead('id', 'FRRRFRRRFWNPNW')
        s = SequenceFeatureAnalysis(read,
                                    landmarkNames=['AlphaHelix', 'BetaTurn'],
                                    trigPointNames=['AminoAcids'])
        self.assertEqual(
            {
                'AlphaHelix': {0, 1, 2, 3, 4, 5, 6, 7, 8},
                'AminoAcids': {9, 13},
                'BetaTurn': {10, 11, 12, 13},
            },
            s.offsets)

    def testPairwiseOffsetOverlaps(self):
        """
        The pairwiseOffsetOverlaps method must return the correct overlap
        matrix.
        """
        read = AARead('id', 'FRRRFRRRFWNPNW')
        s = SequenceFeatureAnalysis(read,
                                    landmarkNames=['AlphaHelix', 'BetaTurn'],
                                    trigPointNames=['AminoAcids'])
        # The order in the overlaps array is AlphaHelix, BetaTurn,
        # AminoAcids, so the 0.2 value corresponds to the
        self.assertTrue(np.array_equal(
            np.array([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.2],
                [0.0, 0.2, 1.0],
            ]),
            s.pairwiseOffsetOverlaps()))

    def testPrintDensities(self):
        """
        The printDensities method must return the correct string.
        """
        read = AARead('id', 'AAFRRRFRRRFWNPNWXX')
        s = SequenceFeatureAnalysis(
            read,
            landmarkNames=['AlphaHelix', 'AlphaHelix_pi', 'BetaTurn'],
            trigPointNames=['AminoAcids'])
        # Note that AlphaHelix_pi does not appear in the output as it finds
        # no features (it will be present in s.landmarksNotFound). Arguably
        # it should also be in the output.
        self.assertEqual(
            (
                'Feature densities:\n'
                '              OVERALL         UNIQUE\n'
                '  AlphaHelix: 50.00% (9/18)   50.00% (9/18)\n'
                '  BetaTurn:   22.22% (4/18)   16.67% (3/18)\n'
                '  AminoAcids: 11.11% (2/18)    5.56% (1/18)\n'
                '\n'
                '22.22% (4/18) of sequence offsets were not covered by any '
                'feature.'
            ),
            s.printDensities())
