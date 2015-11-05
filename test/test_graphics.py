from unittest import TestCase

from light.graphics import PlotHashesInSubjectAndRead

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

    def testOnlyShowBestBin(self):
        """
        The onlyShowBestBin option must work correctly, by setting the 'bins'
        attribute to contain just one bin.
        """

        A = 'FRRRFRRRFXXXXXX'
        C = 'FRRRRFRRRRFXXXXXX'

        subject = AARead('subject', 5 * A + C + 2 * A)
        query = AARead('query', 5 * A)

        significanceFraction = 0.01

        # There are 11 significant bins.
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_pi'],
            trigPointNames=[], distanceBase=1.025, limitPerLandmark=50,
            minDistance=1, maxDistance=100, showInsignificant=False,
            significanceFraction=significanceFraction)
        self.assertEqual(11, len(hashes.bins))

        bestScore = hashes.result.analysis['bestScore']

        # Same input, but restricting ourselves to only the single most
        # significant bin:
        hashes = PlotHashesInSubjectAndRead(
            query, subject, landmarkNames=['AlphaHelix', 'AlphaHelix_pi'],
            trigPointNames=[], distanceBase=1.025, limitPerLandmark=50,
            minDistance=1, maxDistance=100, showInsignificant=False,
            significanceFraction=significanceFraction, onlyShowBestBin=True)
        self.assertEqual(1, len(hashes.bins))

        # Check that the best bin when we use onlyShowBestBin has the same
        # score as the bin we get when we don't use onlyShowBestBin.
        self.assertEqual(bestScore, hashes.result.analysis['bestScore'])
