from unittest import TestCase

from light.landmarks.alpha_helix import AlphaHelix
from light.trig.peaks import Peaks
from light.database import ScannedReadDatabase
from dark.reads import AARead


class TestScannedReadDatabase(TestCase):
    """
    Tests for the light.database.ScannedReadDatabase class.
    """
    def testFindersAreStored(self):
        """
        The list of landmark and trig point finders must be stored correctly.
        """
        db = ScannedReadDatabase([AlphaHelix], [Peaks])
        self.assertEqual([AlphaHelix], db.landmarkFinderClasses)
        self.assertEqual([Peaks], db.trigPointFinderClasses)

    def testInitialStatistics(self):
        """
        The database statistics must be initially correct.
        """
        db = ScannedReadDatabase([], [])
        self.assertEqual(0, db.readCount)
        self.assertEqual(0, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testInitialDatabaseIsEmpty(self):
        """
        The database must be empty if no reads have been added.
        """
        db = ScannedReadDatabase([AlphaHelix], [Peaks])
        self.assertEqual({}, db.d)

    def testOneReadOneLandmark(self):
        """
        If one read is added but it only has one landmark, nothing is added
        to the database.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id', 'FRRRFRRRF'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkStatistics(self):
        """
        If one read is added the database statistics must be correct.
        """
        db = ScannedReadDatabase([], [])
        db.addRead(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.readCount)
        self.assertEqual(9, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testOneReadTwoLandmarks(self):
        """
        If one read is added and it has two landmarks, two keys are added
        to the database.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:-23': set([('id', 0)]),
                'A3:A2:23': set([('id', 23)]),
            },
            db.d)

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one read is added, the database statistics must be correct.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(1, db.readCount)
        self.assertEqual(36, db.totalResidues)
        self.assertEqual(22, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarks(self):
        """
        If two identical reads are added, both with two landmarks, two keys
        are added to the database and both reads are listed in the dictionary
        values for those two keys.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id1', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addRead(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:-23': set([('id1', 0), ('id2', 0)]),
                'A3:A2:23': set([('id2', 23), ('id1', 23)]),
            },
            db.d)

    def testTwoReadsTwoLandmarksStatistics(self):
        """
        If two identical reads are added, the database statistics must be
        correct.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id1', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addRead(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(2, db.readCount)
        self.assertEqual(72, db.totalResidues)
        self.assertEqual(44, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarksDifferentOffsets(self):
        """
        If two reads are added, both with two landmarks separated by the same
        distance, two keys are added to the database and both reads are listed
        in the dictionary values for those two keys.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id1', 'AFRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addRead(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:-23': set([('id1', 1), ('id2', 0)]),
                'A3:A2:23': set([('id2', 23), ('id1', 24)]),
            },
            db.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one read is added and it has one landmark and one peak, one pair is
        added to the database.
        """
        db = ScannedReadDatabase([AlphaHelix], [Peaks])
        db.addRead(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual(
            {
                'A2:P:-10': set([('id', 0)]),
            },
            db.d)

    def testOneReadOneLandmarkOnePeakNoTrigFinders(self):
        """
        If one read is added and it has one landmark and one peak, but no trig
        finders are in use, nothing is added to the database.
        """
        db = ScannedReadDatabase([AlphaHelix], [])
        db.addRead(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual({}, db.d)
