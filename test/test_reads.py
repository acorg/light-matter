from unittest import TestCase

from dark.reads import AARead

from light.reads import ScannedRead
from light.features import Landmark, TrigPoint


class TestScannedRead(TestCase):
    """
    Tests for the light.reads.ScannedRead class
    """
    def testReadIsStored(self):
        """
        A scanned read must store the read it is passed.
        """
        read = AARead('id', 'AAA')
        scannedRead = ScannedRead(read)
        self.assertIs(read, scannedRead.read)

    def testNoLandmarks(self):
        """
        A freshly created scanned read must have no landmarks.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        self.assertEqual(0, len(read.landmarks))

    def testNoTrigPoints(self):
        """
        A freshly created scanned read must have no trig points.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        self.assertEqual(0, len(read.trigPoints))

    def testAddLandmark(self):
        """
        It must be possible to add a landmark to a scanned read.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        landmark = Landmark('name', 'symbol', 0, 1)
        read.landmarks.append(landmark)
        self.assertIs(landmark, read.landmarks[0])

    def testAddTrigPoint(self):
        """
        It must be possible to add a trig point to a scanned read.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        self.assertIs(trigPoint, read.trigPoints[0])

    def testNoCoverage(self):
        """
        If a scanned read has no landmarks or trig points, its coveredIndices
        method must return the empty set.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        self.assertEqual(set(), read.coveredIndices())

    def testCoverageOfOneTrigPoint(self):
        """
        If a scanned read has one trig point, its coveredIndices method must
        return the index of that trig point.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        self.assertEqual({0}, read.coveredIndices())

    def testCoverageOfTwoTrigPoints(self):
        """
        If a scanned read has two trig points, its coveredIndices method must
        return the indices of those trig points.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        trigPoint = TrigPoint('name', 'symbol', 1)
        read.trigPoints.append(trigPoint)
        self.assertEqual({0, 1}, read.coveredIndices())

    def testCoverageOfTwoIdenticalTrigPoints(self):
        """
        If a scanned read has two trig points that occur at the same index, its
        coveredIndices method must return the common index of the trig point.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        self.assertEqual({0}, read.coveredIndices())

    def testCoverageOfOneLandmarkPoint(self):
        """
        If a scanned read has one landmark, its coveredIndices method must
        return the indices of that landmark.
        """
        read = ScannedRead(AARead('id', 'AAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        self.assertEqual({0, 1}, read.coveredIndices())

    def testCoverageOfTwoLandmarkPoints(self):
        """
        If a scanned read has two landmarks, its coveredIndices method must
        return the indices of those two landmarks.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        landmark = Landmark('name', 'symbol', 3, 2)
        read.landmarks.append(landmark)
        self.assertEqual({0, 1, 3, 4}, read.coveredIndices())
