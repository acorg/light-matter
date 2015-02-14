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

    def testCoverageOfTwoDifferentNotOverlappingLandmarkPoints(self):
        """
        If a scanned read has two landmarks, its coveredIndices method must
        return the indices of those two landmarks.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol1', 0, 2)
        read.landmarks.append(landmark)
        landmark = Landmark('name', 'symbol2', 4, 2)
        read.landmarks.append(landmark)
        self.assertEqual({0, 1, 4, 5}, read.coveredIndices())

    def testCoverageOfTwoDifferentOverlappingLandmarkPoints(self):
        """
        If a scanned read has two landmarks, its coveredIndices method must
        return the indices of those two landmarks.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol1', 0, 5)
        read.landmarks.append(landmark)
        landmark = Landmark('name', 'symbol2', 4, 2)
        read.landmarks.append(landmark)
        self.assertEqual({0, 1, 2, 3, 4, 5}, read.coveredIndices())

    def testGetPairsNoLandmarksNoTrigPoints(self):
        """
        If a scanned read has no landmarks or trig points, its getPairs
        method should generate no pairs.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        self.assertEqual([], list(read.getPairs()))

    def testGetPairsOneLandmark(self):
        """
        If a scanned read has one landmark, its getPairs method should generate
        no pairs because there is nothing to pair with.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        self.assertEqual([], list(read.getPairs()))

    def testGetPairsOneLandmarkOneTrigPoint(self):
        """
        If a scanned read has one landmark and one trig point, its getPairs
        method should generate one pair.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        self.assertEqual([(landmark, trigPoint)], list(read.getPairs()))

    def testGetPairsOneLandmarkOneTrigPointLimitZero(self):
        """
        If a scanned read has one landmark and one trig point, but getPairs
        is given a limit of 0, no pairs should be generated.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        trigPoint = TrigPoint('name', 'symbol', 0)
        read.trigPoints.append(trigPoint)
        self.assertEqual([], list(read.getPairs(limitPerLandmark=0)))

    def testGetPairsOneLandmarkTwoTrigPoints(self):
        """
        If a scanned read has one landmark and two trig points, its getPairs
        method should generate two pairs.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark = Landmark('name', 'symbol', 0, 2)
        read.landmarks.append(landmark)
        trigPoint1 = TrigPoint('name', 'symbol', 1)
        trigPoint2 = TrigPoint('name', 'symbol', 2)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs())
        self.assertEqual([(landmark, trigPoint1), (landmark, trigPoint2)],
                         result)

    def testGetPairsTwoLandmarksOneTrigPoint(self):
        """
        If a scanned read has two landmarks and one trig point, its getPairs
        method should generate four pairs in the correct order of closeness,
        not yielding (landmark, landmark) pairs with the second landmark to
        the left of the first until all other pairs have been yielded.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'symbol', 0, 2)
        landmark2 = Landmark('name', 'symbol', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint = TrigPoint('name', 'symbol', 4)
        read.trigPoints.append(trigPoint)
        result = list(read.getPairs())
        self.assertEqual(
            [(landmark1, landmark2), (landmark1, trigPoint),
             (landmark2, trigPoint), (landmark2, landmark1)],
            result)

    def testGetPairsTwoLandmarksTwoTrigPoints(self):
        """
        If a scanned read has two landmarks and two trig points, its getPairs
        method should generate six pairs, in the correct order of closeness,
        not yielding (landmark, landmark) pairs with the second landmark to
        the left of the first until all other pairs have been yielded.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'symbol', 0, 2)
        landmark2 = Landmark('name', 'symbol', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint1 = TrigPoint('name', 'symbol', 4)
        trigPoint2 = TrigPoint('name', 'symbol', 5)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs())
        self.assertEqual(
            [(landmark1, landmark2), (landmark1, trigPoint1),
             (landmark1, trigPoint2),
             (landmark2, trigPoint1), (landmark2, trigPoint2),
             (landmark2, landmark1)],
            result)

    def testGetPairsTwoLandmarksTwoTrigPointsLimitOne(self):
        """
        If a scanned read has two landmarks and two trig points, its getPairs
        method should generate two pairs if a limit of 1 nearby feature per
        landmark is passed.

        The (landmark2, landmark1) pair is not produced due to the
        limitPerLandmark=1 value. This is because the getPairs method gives
        preference to features that occur to the right of landmarks. In this
        case, even though landmark1 is closer to landmark2 than trigPoint1 is,
        the (landmark2, trigPoint1) pair is yielded because trigPoint1 is to
        the right of landmark2. The (landmark2, trigPoint2) and finally the
        (landmark2, landmark1) pairs would be the next to be yielded, but the
        limitPerLandmark limit has already been reached.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'L1', 0, 2)
        landmark2 = Landmark('name', 'L2', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint1 = TrigPoint('name', 'T1', 4)
        trigPoint2 = TrigPoint('name', 'T2', 5)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs(limitPerLandmark=1))
        self.assertEqual([(landmark1, landmark2), (landmark2, trigPoint1)],
                         result)

    def testGetPairsTwoLandmarksTwoTrigPointsLimitTwo(self):
        """
        If a scanned read has two landmarks and two trig points, its getPairs
        method should generate four pairs if a limit of 2 nearby features per
        landmark is passed.

        The (landmark2, landmark1) pair is not produced due to the
        limitPerLandmark=2 value. This is because the getPairs method gives
        preference to features that occur to the right of landmarks. In this
        case, even though landmark1 is closer to landmark2 than trigPoint2 is,
        the (landmark2, trigPoint2) pair is yielded because trigPoint2 is to
        the right of landmark2. The (landmark2, landmark1) pair would be the
        next to be yielded, but by then the limitPerLandmark limit has been
        reached.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'L1', 0, 2)
        landmark2 = Landmark('name', 'L2', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint1 = TrigPoint('name', 'T1', 4)
        trigPoint2 = TrigPoint('name', 'T2', 5)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs(limitPerLandmark=2))
        self.assertEqual(
            [(landmark1, landmark2), (landmark1, trigPoint1),
             (landmark2, trigPoint1), (landmark2, trigPoint2)],
            result)

    def testGetPairsTwoLandmarksTwoTrigPointsMaxDistanceThree(self):
        """
        If a scanned read has two landmarks and two trig points, its getPairs
        method should generate four pairs if a maximum distance of 3 is
        passed and one trig point is too far away.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'L1', 0, 2)
        landmark2 = Landmark('name', 'L2', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint1 = TrigPoint('name', 'T1', 3)
        trigPoint2 = TrigPoint('name', 'T2', 5)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs(maxDistance=3))
        self.assertEqual(
            [(landmark1, landmark2), (landmark1, trigPoint1),
             (landmark2, trigPoint1), (landmark2, landmark1)],
            result)

    def testGetPairsTwoLandmarksTwoTrigPointsMinDistanceThree(self):
        """
        If a scanned read has two landmarks and two trig points, its getPairs
        method should generate one pair if a minimum distance of 3 is
        passed and three trig points are too close.
        """
        read = ScannedRead(AARead('id', 'AAAAA'))
        landmark1 = Landmark('name', 'L1', 0, 2)
        landmark2 = Landmark('name', 'L2', 1, 2)
        read.landmarks.extend([landmark1, landmark2])
        trigPoint1 = TrigPoint('name', 'T1', 3)
        trigPoint2 = TrigPoint('name', 'T2', 5)
        read.trigPoints.extend([trigPoint1, trigPoint2])
        result = list(read.getPairs(minDistance=5))
        self.assertEqual([(landmark1, trigPoint2)], result)
