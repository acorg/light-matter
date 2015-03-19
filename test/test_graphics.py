from unittest import TestCase

from light.graphics import PlotHashesInSubjectAndRead

from dark.reads import AARead


class TestPlotHashesInSubjectAndRead(TestCase):
    """
    Tests for the light.graphics.PlotHashesInSubjectAndRead class.
    """
    def testNoHashes(self):
        """
        If there are no hashes in subject and query, matchingFeatures,
        queryFeatures and subjectFeatures must be empty.
        """
        subject = AARead('subject', 'RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR')
        query = AARead('query', 'AAAAAAAAAAAAAAAAAAAAAAA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.25,
                                            maxDistance=10, minDistance=1,
                                            limitPerLandmark=10,
                                            distanceBase=1.0)
        self.assertEqual(0, len(list(hashes.matchingFeatures)))
        self.assertEqual(0, len(list(hashes.queryFeatures)))
        self.assertEqual(0, len(list(hashes.subjectFeatures)))

    def testNoMatchingHashes(self):
        """
        If there are no matching hashes in subject and query, matchingFeatures
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
                                            distanceBase=1.0)
        self.assertEqual(0, len(list(hashes.matchingFeatures)))
        self.assertEqual(1, len(list(hashes.queryFeatures)))
        self.assertEqual(1, len(list(hashes.subjectFeatures)))

    def testNoQuerysubjectHashes(self):
        """
        If all subjectFeatures in the query also occur in the subject,
        queryFeatures must be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAASARRRRFRRFRRFAAASA')
        query = AARead('query', 'ASARRRRFRRFRRFAAASA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix',
                                                           'AlphaHelix_3_10'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.01,
                                            maxDistance=50, minDistance=1,
                                            limitPerLandmark=10,
                                            distanceBase=1.0)
        matchingHashesCount = sum([len(bin_) for bin_ in
                                   hashes.matchingFeatures])
        self.assertEqual(2, matchingHashesCount)
        self.assertEqual(0, len(hashes.queryFeatures))
        self.assertEqual(3, len(hashes.subjectFeatures))

    def testNoSubjectHashes(self):
        """
        If all hashes in the subject also occur in the query, subjectFeatures
        must be empty.
        """
        subject = AARead('subject', 'FRRFRRFRRFAAAAAAAAAAAASA')
        query = AARead('query', 'FRRFRRFRRFAAAAAAAAAAAASARRRRFRRFRRFAAASA')
        hashes = PlotHashesInSubjectAndRead(query, subject,
                                            landmarkNames=['AlphaHelix',
                                                           'AlphaHelix_3_10'],
                                            trigPointNames=['Peaks'],
                                            significanceFraction=0.01,
                                            maxDistance=50, minDistance=1,
                                            limitPerLandmark=10,
                                            distanceBase=1.0)
        matchingHashesCount = sum([len(bin_) for bin_ in
                                   hashes.matchingFeatures])
        self.assertEqual(1, matchingHashesCount)
        self.assertEqual(4, len(hashes.queryFeatures))
        self.assertEqual(0, len(hashes.subjectFeatures))
