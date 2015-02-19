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
        queryHashes and subjectHashes must be empty.
        """
        subject = AARead('subject', 'RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR')
        query = AARead('query', 'AAAAAAAAAAAAAAAAAAAAAAA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.25,
                                            maxDistance=10, minDistance=1,
                                            limitPerLandmark=10,
                                            bucketFactor=1)
        self.assertEqual(0, len(hashes.matchingHashes))
        self.assertEqual(0, len(hashes.queryHashes))
        self.assertEqual(0, len(hashes.subjectHashes))

    def testNoMatchingHashes(self):
        """
        If there are no matching hashes in subject and query, matchingHashes
        must be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAAASARRRRRRRRRRRRRR')
        query = AARead('query', 'FRRFRRFAAASAAAAAAAAAAAAA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix',
                                                           'AlphaHelix_3_10'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.25,
                                            maxDistance=50, minDistance=1,
                                            limitPerLandmark=10,
                                            bucketFactor=1)
        self.assertEqual(0, len(hashes.matchingHashes))
        self.assertEqual(1, len(hashes.queryHashes))
        self.assertEqual(1, len(hashes.subjectHashes))

    def testNoQueryHashes(self):
        """
        If all hashes in the query also occur in the subject, queryHashes must
        be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAASARRRRFRRFRRFAAASA')
        query = AARead('query', 'ASARRRRFRRFRRFAAASA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix',
                                                           'AlphaHelix_3_10'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.25,
                                            maxDistance=50, minDistance=1,
                                            limitPerLandmark=10,
                                            bucketFactor=1)
        self.assertEqual(0, len(hashes.matchingHashes))
        self.assertEqual(1, len(hashes.queryHashes))
        self.assertEqual(1, len(hashes.subjectHashes))

    def testNoSubjectHashes(self):
        """
        If all hashes in the subject also occur in the query, subjectHashes
        must be empty.
        """
        pass
