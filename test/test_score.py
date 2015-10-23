import warnings
from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark, TrigPoint
from light.parameters import Parameters
from light.subject import Subject
from light.score import (
    MinHashesScore, FeatureMatchingScore, histogramBinFeatures)
from light.histogram import Histogram
from light.landmarks import AlphaHelix


class TestMinHashesScore(TestCase):
    """
    Tests for the light.score.MinHashesScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0
        """
        histogram = Histogram()
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        self.assertEqual(0.0, mhs.calculateScore(0))

    def testZeroMinHashCount(self):
        """
        If the minimum hash count for the query/subject is zero, the score
        for any bin must be 0.0
        """
        # Use a histogram with one bin, and put something in it, to make sure
        # we're not testing an empty histogram bin.
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        self.assertEqual(0.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        one thing in it, the score for the bin must be 1.0
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        self.assertEqual(1.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountTwo(self):
        """
        If the minimum hash count for the query/subject is two and the bin has
        one thing in it, the score for the bin must be 0.5
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 2)
        self.assertEqual(0.5, mhs.calculateScore(0))

    def testBinWithTwoElementsMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        two things in it, a warning will be issued to say that the bin has more
        things in it than should be possible (given the minimum hash count of
        the query/subject).
        """
        histogram = Histogram(1)
        histogram.add(444)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            mhs.calculateScore(0)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 2 deltas for a query/subject pair with a '
                     'minimum hash count of only 1.')
            self.assertIn(error, str(w[0].message))


class TestHistogramBinFeatures(TestCase):
    """
    Tests for the light.score.histogramBinFeatures function.
    """
    def testInvalidQueryOrSubjectSpecifier(self):
        """
        If an invalid value for queryOrSubject ('xxx') is passed to
        histogramBinFeatures, it must raise a KeyError (when looking
        in the bin item dictionary for 'xxxOffsets').
        """
        histogram = Histogram(1)
        histogram.add(44, {
            # It doesn't matter what values we pass here, as the KeyError
            # will be triggered by the lack of a dictionary key.
            'landmark': None,
            'queryOffsets': None,
            'subjectOffsets': None,
            'trigPoint': None,
        })
        histogram.finalize()
        self.assertRaisesRegexp(KeyError, 'xxxOffsets', histogramBinFeatures,
                                histogram[0], 'xxx')

    def testEmptyBin(self):
        """
        If a histogram bin is empty, then when asked to extract its features
        and offsets, histogramBinFeatures must return two empty sets.
        """
        histogram = Histogram()
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        self.assertEqual(set(), features)
        self.assertEqual(set(), offsets)

    def testOneFeatureQuery(self):
        """
        If a histogram bin has just one feature, histogramBinFeatures must
        return the details of that feature in the query.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110]],
            'subjectOffsets': None,
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        self.assertEqual(set([landmark, trigPoint]), features)
        self.assertEqual(set([100, 110, 120]), offsets)

    def testOneFeatureSubject(self):
        """
        If a histogram bin has just one feature, histogramBinFeatures must
        return the details of that feature in the subject.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': None,
            'subjectOffsets': [[100, 110]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'subject')
        self.assertEqual(set([landmark, trigPoint]), features)
        self.assertEqual(set([100, 110, 120]), offsets)

    def testOneFeatureTwoLocations(self):
        """
        If a histogram bin has one feature that appears in two places,
        histogramBinFeatures must return the details of that feature.
        """
        landmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint1 = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark1,
            'queryOffsets': [[100, 110], [200, 210]],
            'subjectOffsets': None,
            'trigPoint': trigPoint1,
        })
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        # Because the hash occurs at two locations, we get back 2 landmarks
        # and 2 trig points from histogramBinFeatures.
        landmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        trigPoint2 = TrigPoint('Peak', 'P', 210)
        self.assertEqual(set([landmark1, trigPoint1, landmark2, trigPoint2]),
                         features)
        self.assertEqual(set([100, 110, 120, 200, 210, 220]), offsets)


class TestFeatureMatchingScore(TestCase):
    """
    Tests for the light.score.FeatureMatchingScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0 if the query and
        subject both have no features.
        """
        histogram = Histogram()
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(0.0, fms.calculateScore(0))

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four) if the query and subject have no
        additional (non-matching) features. The reward is multiplied by four
        because both the landmark and the trig point (2) appear in both the
        query and the subject (2) and 2 x 2 = 4.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110]],
            'subjectOffsets': [[100, 110]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinOccurringInTwoPlaces(self):
        """
        A bin containing one hash that occurs in two place must have a score
        that is the feature match reward multiplied by eight) if the query and
        subject have no additional (non-matching) features. The reward is
        multiplied by eight because both the landmark and the trig point (2)
        appear twice (2) in both the query and the subject (2) and
        2 x 2 x 2 = 8.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110], [200, 210]],
            'subjectOffsets': [[100, 110], [200, 210]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(8 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasOneFeatureOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has an additional feature
        that is outside the match area.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110]],
            'subjectOffsets': [[100, 110]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasTwoFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has two additional
        features that are outside the match area.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110]],
            'subjectOffsets': [[100, 110]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF' + ('F' * 200) + 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryAndSubjectHaveOneFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query and the subject each have
        one additional feature that are both outside the match area.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        trigPoint = TrigPoint('Peak', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[100, 110]],
            'subjectOffsets': [[100, 110]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', ('F' * 200) + 'FRRRFRRRF', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area.
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 0, 20)
        trigPoint = TrigPoint('Peak', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[0, 10]],
            'subjectOffsets': [[0, 10]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE +
            Parameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExactlySpanningMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that exactly spans the match area.
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 0, 9)
        trigPoint = TrigPoint('Peak', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[0, 5]],
            'subjectOffsets': [[0, 5]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE +
            Parameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExceedingMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        exceeds the match area on both sides.
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 2, 5)
        trigPoint = TrigPoint('Peak', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[2, 5]],
            'subjectOffsets': [[2, 5]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasTwoUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus twice the feature mismatch score if
        the query has an additional feature that is inside the match area.
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 0, 20)
        trigPoint = TrigPoint('Peak', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[0, 10]],
            'subjectOffsets': [[0, 10]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRFAFRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE +
            2 * Parameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the left).
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 5, 20)
        trigPoint = TrigPoint('Peak', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[5, 10]],
            'subjectOffsets': [[5, 10]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the right).
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 2, 3)
        trigPoint = TrigPoint('Peak', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[2, 5]],
            'subjectOffsets': [[2, 5]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'AAAFRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(
            4 * Parameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureInMatchNonDefault(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area,
        including when non-default values are used for the feature match and
        mismatch scores.
        """
        landmark = Landmark('AlphaHelix_pi', 'C', 0, 20)
        trigPoint = TrigPoint('Peak', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'landmark': landmark,
            'queryOffsets': [[0, 10]],
            'subjectOffsets': [[0, 10]],
            'trigPoint': trigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [], featureMatchScore=3.1,
                            featureMismatchScore=-1.2)
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        self.assertEqual(4 * 3.1 - 1.2, fms.calculateScore(0))
