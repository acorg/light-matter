from unittest import TestCase
from random import shuffle as pyshuffle
from copy import copy

from light.features import (
    _Feature, Landmark, TrigPoint, CombinedFeatureList, Finder)


def shuffle(l):
    """
    Copy a list, shuffle it, and return it.

    @param l: A C{list}.
    @return: A new list, with the same elements, but in a randomly shuffled
        order.
    """
    result = copy(l)
    pyshuffle(result)
    return result


class TestFeature(TestCase):
    """
    Tests for the light.features._Feature class
    """
    def testZeroCoverage(self):
        """
        A feature of zero length must not cover any AA offsets.
        """
        self.assertEqual(set(),
                         _Feature('name', 'symbol', 10, 0).coveredOffsets())

    def testCoverage(self):
        """
        The coveredOffsets function must give the expected result for a
        non-zero length feature.
        """
        self.assertEqual(set([10, 11, 12, 13, 14]),
                         _Feature('name', 'symbol', 10, 5).coveredOffsets())


class TestLandmarks(TestCase):
    """
    Tests for the light.features.Landmark class
    """
    def testEqual(self):
        """
        Identical landmarks must compare equal.
        """
        landmark1 = Landmark('name', 'L', 0, 1)
        landmark2 = Landmark('name', 'L', 0, 1)
        self.assertEqual(landmark1, landmark2)

    def testDifferingNamesNonEqual(self):
        """
        Landmarks with different names must not compare equal.
        """
        landmark1 = Landmark('name1', 'L', 0, 1)
        landmark2 = Landmark('name2', 'L', 0, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingSymbolsNonEqual(self):
        """
        Landmarks with different symbols must not compare equal.
        """
        landmark1 = Landmark('name', 'L1', 0, 1)
        landmark2 = Landmark('name', 'L2', 0, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingOffsetsNonEqual(self):
        """
        Landmarks with different offsets must not compare equal.
        """
        landmark1 = Landmark('name', 'L', 0, 1)
        landmark2 = Landmark('name', 'L', 1, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingLengthsNonEqual(self):
        """
        Landmarks with different lengths must not compare equal.
        """
        landmark1 = Landmark('name', 'L', 0, 1)
        landmark2 = Landmark('name', 'L', 0, 2)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingSymbolDetailsNonEqual(self):
        """
        Landmarks with different symbol details must not compare equal.
        """
        landmark1 = Landmark('name', 'L', 0, 1, 0)
        landmark2 = Landmark('name', 'L', 0, 1, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testHashkeyNoSymbolDetail(self):
        """
        The hashkey function must return as expected when no value is given for
        the symbol detail.
        """
        landmark = Landmark('name', 'L', 0, 1)
        self.assertEqual('L', landmark.hashkey())

    def testHashkeyWithSymbolDetail(self):
        """
        The hashkey function must return as expected when symbol detail is
        given.
        """
        landmark = Landmark('name', 'L', 0, 1, 2)
        self.assertEqual('L2', landmark.hashkey())

    def testSortOnName(self):
        """
        Sorting must be according to name if all else is equal.
        """
        expected = [
            Landmark('nameA', 'L', 50, 22, 0),
            Landmark('nameB', 'L', 50, 22, 0),
            Landmark('nameC', 'L', 50, 22, 0),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testSortOnSymbol(self):
        """
        Sorting does not consider the symbol becaused it uses the name.
        All names are unique and each has a unique symbol, so nothing
        would be gained by additionally sorting on symbol once sorting
        on name had been done.
        """
        pass

    def testSortOnOffset(self):
        """
        Sorting must be according to offset if all else is equal.
        """
        expected = [
            Landmark('name', 'L', 50, 22, 0),
            Landmark('name', 'L', 51, 22, 0),
            Landmark('name', 'L', 52, 22, 0),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testSortOnLength(self):
        """
        Sorting must be according to length if all else is equal.
        """
        expected = [
            Landmark('name', 'L', 50, 20, 0),
            Landmark('name', 'L', 50, 21, 0),
            Landmark('name', 'L', 50, 22, 0),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testSortOnSymbolDetail(self):
        """
        Sorting must be according to symbol detail if all else is equal.
        """
        expected = [
            Landmark('name', 'L', 50, 22, 0),
            Landmark('name', 'L', 51, 22, 1),
            Landmark('name', 'L', 52, 22, 2),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testLandmarkCoveredFeatureIndices(self):
        """
        The right offsets covered by this landmark must be returned.
        """
        landmark = Landmark('name', 'L', 50, 5, 0)
        self.assertEqual({50, 51, 52, 53, 54},
                         landmark.coveredFeatureIndices())


class TestTrigPoints(TestCase):
    """
    Tests for the light.features.TrigPoint class
    """
    def testEqual(self):
        """
        Identical trig points must compare equal.
        """
        trigPoint1 = TrigPoint('name', 'T', 0)
        trigPoint2 = TrigPoint('name', 'T', 0)
        self.assertEqual(trigPoint1, trigPoint2)

    def testDifferingNamesNonEqual(self):
        """
        Trig points with different names must not compare equal.
        """
        trigPoint1 = TrigPoint('name1', 'T', 0)
        trigPoint2 = TrigPoint('name2', 'T', 0)
        self.assertNotEqual(trigPoint1, trigPoint2)

    def testDifferingSymbolsNonEqual(self):
        """
        Trig points with different symbols must not compare equal.
        """
        trigPoint1 = TrigPoint('name', 'T1', 0)
        trigPoint2 = TrigPoint('name', 'T2', 0)
        self.assertNotEqual(trigPoint1, trigPoint2)

    def testDifferingOffsetsNonEqual(self):
        """
        Trig points with different offsets must not compare equal.
        """
        trigPoint1 = TrigPoint('name', 'T', 0)
        trigPoint2 = TrigPoint('name', 'T', 1)
        self.assertNotEqual(trigPoint1, trigPoint2)

    def testTrigPointHashkey(self):
        """
        The hashkey function must return as expected.
        """
        trigPoint = TrigPoint('name', 'L', 0)
        self.assertEqual('L', trigPoint.hashkey())

    def testSortOnName(self):
        """
        Sorting must be according to name if all else is equal.
        """
        expected = [
            TrigPoint('nameA', 'L', 50),
            TrigPoint('nameB', 'L', 50),
            TrigPoint('nameC', 'L', 50),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testSortOnSymbol(self):
        """
        Sorting does not consider the symbol becaused it uses the name.
        All names are unique and each has a unique symbol, so nothing
        would be gained by additionally sorting on symbol once sorting
        on name had been done.
        """
        pass

    def testSortOnOffset(self):
        """
        Sorting must be according to offset if all else is equal.
        """
        expected = [
            TrigPoint('name', 'L', 50),
            TrigPoint('name', 'L', 51),
            TrigPoint('name', 'L', 52),
        ]

        for _ in range(100):
            self.assertEqual(expected, sorted(shuffle(expected)))

    def testTrigPointCoveredFeatureIndices(self):
        """
        The right offset covered by this trigPoint must be returned.
        """
        trigPoint = TrigPoint('name', 'L', 50)
        self.assertEqual({50}, trigPoint.coveredFeatureIndices())


class TestCombinedFeatureList(TestCase):
    """
    Tests for the light.features.CombinedFeatureList class
    """
    def testEmpty(self):
        """
        If a combined feature list has no landmarks or trig points, its nearest
        method should not generate any features.
        """
        cfl = CombinedFeatureList([], [])
        self.assertEqual([], list(cfl.nearest(0)))

    def testNothingInRange(self):
        """
        If a combined feature list has no landmarks or trig points within range
        of the offset passed to its nearest method, the method should not
        generate any features.
        """
        landmark1 = Landmark('name1', 'L1', 0, 1)
        landmark2 = Landmark('name2', 'L2', 10, 1)
        cfl = CombinedFeatureList([landmark1, landmark2], [])
        self.assertEqual([], list(cfl.nearest(5, maxDistance=2)))

    def testNearestFindsZeroDistanceLandmark(self):
        """
        If a combined feature list has one landmark and the offset of that
        landmark is given to its nearest method, that landmark should be
        generated.
        """
        landmark = Landmark('name', 'L', 10, 1)
        cfl = CombinedFeatureList([landmark], [])
        result = list(cfl.nearest(10, maxDistance=0))
        self.assertEqual(1, len(result))
        self.assertIs(landmark, result[0])

    def testNearestFindsZeroDistanceTrigPoint(self):
        """
        If a combined feature list has one trig point and the offset of that
        trig point is given to its nearest method, that trig point should be
        generated.
        """
        trigPoint = TrigPoint('name', 'T', 10)
        cfl = CombinedFeatureList([], [trigPoint])
        result = list(cfl.nearest(10, maxDistance=0))
        self.assertEqual(1, len(result))
        self.assertIs(trigPoint, result[0])

    def testFindLandmarkAndTrigPoint(self):
        """
        If a combined feature list has one landmark and one trig point within
        range, with the landmark closest to the passed offset, the nearest
        method should yield them both, with the landmark first.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(4))
        self.assertEqual([landmark, trigPoint], result)

    def testFindTrigPointAndLandmarkWithLandmarkFirst(self):
        """
        If a combined feature list has one landmark and one trig point within
        range, with the trig point closest to the passed offset, the nearest
        method should yield them both, with the landmark first (because
        landmarks are prioritized ahead of trig points).
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(6))
        self.assertEqual([landmark, trigPoint], result)

    def testFindLandmarkAndTrigPointMaxDistance(self):
        """
        If a combined feature list has one landmark and one trig point
        with only the landmark within the passed max distance, the nearest
        method should yield just the landmark.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(4, maxDistance=5))
        self.assertEqual([landmark], result)

    def testFindTrigPointAndLandmarkMaxDistance(self):
        """
        If a combined feature list has one landmark and one trig point
        with only the trig point within the passed max distance, the nearest
        method should yield just the trig point.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(6, maxDistance=5))
        self.assertEqual([trigPoint], result)

    def testFindLandmarkAndTrigPointWithOffsetBiggerThanBoth(self):
        """
        If a combined feature list has one landmark and one trig point and the
        passed offset is larger than both their offsets, the nearest method
        should yield them both (landmarks first, due to landmark priority,
        even though it is further away than the trig point).
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(20))
        self.assertEqual([landmark, trigPoint], result)

    def testFindMany(self):
        """
        If a combined feature list has many landmarks and trig points
        the nearest method should yield them all.
        """
        landmark1 = Landmark('name1', 'L1', 0, 1)
        landmark2 = Landmark('name2', 'L2', 3, 1)
        landmark3 = Landmark('name3', 'L3', 5, 1)
        landmark4 = Landmark('name4', 'L4', 7, 1)
        trigPoint1 = TrigPoint('name5', 'T1', 9)
        trigPoint2 = TrigPoint('name6', 'T2', 11)
        trigPoint3 = TrigPoint('name7', 'T3', 13)
        cfl = CombinedFeatureList(
            [landmark1, landmark2, landmark3, landmark4],
            [trigPoint1, trigPoint2, trigPoint3])
        result = list(cfl.nearest(0))
        self.assertEqual(
            [landmark1, landmark2, landmark3, landmark4,
             trigPoint1, trigPoint2, trigPoint3],
            result)

    def testFindAllMinDistance(self):
        """
        If a combined feature list has all landmarks and trig points beyond the
        minDistance, find all, with the landmarks first.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint1 = TrigPoint('name2', 'T', 10)
        trigPoint2 = TrigPoint('name3', 'T', 11)
        cfl = CombinedFeatureList([landmark], [trigPoint1, trigPoint2])
        result = list(cfl.nearest(6, minDistance=1))
        self.assertEqual([landmark, trigPoint1, trigPoint2], result)

    def testFindNoneMinDistance(self):
        """
        If a combined feature list has no landmarks and trig points beyond the
        minDistance, find none.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(6, minDistance=10))
        self.assertEqual([], result)

    def testFindOneMinDistance(self):
        """
        If a combined feature list has one of two landmarks and trig points
        inside the minDistance, find one.
        """
        landmark = Landmark('name1', 'L', 0, 1)
        trigPoint = TrigPoint('name2', 'T', 10)
        cfl = CombinedFeatureList([landmark], [trigPoint])
        result = list(cfl.nearest(8, minDistance=5))
        self.assertEqual([landmark], result)

    def testFindOneMinDistanceMaxDistance(self):
        """
        If a combined feature list has three landmarks and trig points, one
        outside the minDistance, one outside the maxDistance, find one.
        """
        trigPoint1 = TrigPoint('name2', 'T', 0)
        landmark = Landmark('name1', 'L', 5, 1)
        trigPoint2 = TrigPoint('name3', 'T', 8)
        cfl = CombinedFeatureList([landmark], [trigPoint1, trigPoint2])
        result = list(cfl.nearest(5, maxDistance=4, minDistance=1))
        self.assertEqual([trigPoint2], result)

    def testFindManyMinDistanceMaxDistance(self):
        """
        If a combined feature list has many landmarks and trig points, some
        outside the minDistance, some outside the maxDistance, find the right
        ones.
        """
        landmark1 = Landmark('name1', 'L1', 0, 1)
        landmark2 = Landmark('name2', 'L2', 3, 1)
        landmark3 = Landmark('name3', 'L3', 5, 1)
        landmark4 = Landmark('name4', 'L4', 7, 1)
        trigPoint1 = TrigPoint('name5', 'T1', 9)
        trigPoint2 = TrigPoint('name6', 'T2', 11)
        trigPoint3 = TrigPoint('name7', 'T3', 13)
        cfl = CombinedFeatureList(
            [landmark1, landmark2, landmark3, landmark4],
            [trigPoint1, trigPoint2, trigPoint3])
        result = list(cfl.nearest(0, maxDistance=10, minDistance=4))
        self.assertEqual(
            [landmark3, landmark4, trigPoint1], result)


class TestSets(TestCase):
    """
    Test set operations on the light.features.Landmark and
    light.features.trigPoint classes.
    """
    def testIdentity(self):
        """
        A set with just one landmark in it should know that it contains that
        landmark.
        """
        landmark = Landmark('name', 'n', 33, 100)
        s = set([landmark])
        self.assertIn(landmark, s)

    def testTuple(self):
        """
        A set with a (landmark, trigPoint) tuple in it should know that it
        contains that tuple.
        """
        landmark = Landmark('name', 'n', 33, 100)
        trigPoint = TrigPoint('trig', 't', 44)
        s = set([(landmark, trigPoint)])
        self.assertIn((landmark, trigPoint), s)


class TestFinder(TestCase):
    """
    Tests for the light.features.Finder class.
    """
    def testDefaultFeatureLengthBase(self):
        """
        An instance must have the right default feature length base.
        """
        finder = Finder()
        self.assertEqual(1.35, finder._featureLengthBase)

    def testNonDefaultFeatureLengthBase(self):
        """
        An instance made with a non-default feature length base must have the
        right featureLengthBase set.
        """
        finder = Finder(featureLengthBase=1.5)
        self.assertEqual(1.5, finder._featureLengthBase)
