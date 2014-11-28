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
