import argparse
from unittest import TestCase
from json import loads
from io import StringIO

from light.parameters import Parameters, FindParameters
from light.landmarks import DEFAULT_LANDMARK_CLASSES, AlphaHelix, BetaStrand
from light.trig import DEFAULT_TRIG_CLASSES, Peaks, Troughs


class TestFindParameters(TestCase):
    """
    Tests for the light.database.FindParameters class.
    """
    def testDefaults(self):
        """
        If no specific parameter values are given, the defaults must be set.
        """
        findParams = FindParameters([], [])
        self.assertEqual(FindParameters.DEFAULT_SCORE_METHOD,
                         findParams.scoreMethod)

    def testNotDefaults(self):
        """
        If specific parameter values are given, the passed values must be set.
        """
        findParams = FindParameters(
            significanceMethod='yyy', significanceFraction=0.5,
            scoreMethod='xxx', featureMatchScore=3.4, featureMismatchScore=9.3)
        self.assertEqual('yyy', findParams.significanceMethod)
        self.assertEqual(0.5, findParams.significanceFraction)
        self.assertEqual('xxx', findParams.scoreMethod)
        self.assertEqual(3.4, findParams.featureMatchScore)
        self.assertEqual(9.3, findParams.featureMismatchScore)

    def testArgs(self):
        """
        It must be possible to parse command line arguments to create a new
        instance of FindParameters.
        """
        parser = argparse.ArgumentParser()
        FindParameters.addArgsToParser(parser)
        args = parser.parse_args([
            '--significanceMethod', 'Always',
            '--significanceFraction', '0.4',
            '--scoreMethod', 'MinHashesScore',
            '--featureMatchScore', '5',
            '--featureMismatchScore', '6',
            '--weights', 'AlphaHelix 2',
        ])

        # Parsing must do the expected thing.
        self.assertEqual('Always', args.significanceMethod)
        self.assertEqual(0.4, args.significanceFraction)
        self.assertEqual('MinHashesScore', args.scoreMethod)
        self.assertEqual(5, args.featureMatchScore)
        self.assertEqual(6, args.featureMismatchScore)

        # We must be able to make an instance from the parsed args.
        findParams = FindParameters.fromArgs(args)
        self.assertEqual('Always', findParams.significanceMethod)
        self.assertEqual(0.4, findParams.significanceFraction)
        self.assertEqual('MinHashesScore', findParams.scoreMethod)
        self.assertEqual(5, findParams.featureMatchScore)
        self.assertEqual(6, findParams.featureMismatchScore)


class TestParameters(TestCase):
    """
    Tests for the light.database.Parameters class.
    """
    def testDefaults(self):
        """
        If no specific parameter values are given, the defaults must be set.
        """
        params = Parameters([], [])
        self.assertEqual([], params.landmarkClasses)
        self.assertEqual([], params.trigPointClasses)
        self.assertEqual(Parameters.DEFAULT_DISTANCE_BASE, params.distanceBase)
        self.assertEqual(Parameters.DEFAULT_LIMIT_PER_LANDMARK,
                         params.limitPerLandmark)
        self.assertEqual(Parameters.DEFAULT_MAX_DISTANCE, params.maxDistance)
        self.assertEqual(Parameters.DEFAULT_MIN_DISTANCE, params.minDistance)

    def testNotDefaults(self):
        """
        If specific parameter values are given, the passed values must be set.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            distanceBase=3.0, limitPerLandmark=10,
                            maxDistance=77, minDistance=66)
        self.assertEqual([AlphaHelix, BetaStrand], params.landmarkClasses)
        self.assertEqual([Peaks, Troughs], params.trigPointClasses)
        self.assertEqual(3.0, params.distanceBase)
        self.assertEqual(10, params.limitPerLandmark)
        self.assertEqual(77, params.maxDistance)
        self.assertEqual(66, params.minDistance)

    def testDefaultLandmarkClasses(self):
        """
        If None is passed for landmark classes, the default list of classes
        must be used.
        """
        params = Parameters(None, [])
        self.assertEqual(DEFAULT_LANDMARK_CLASSES, params.landmarkClasses)

    def testDefaultTrigPointClasses(self):
        """
        If None is passed for trig point classes, the default list of classes
        must be used.
        """
        params = Parameters([], None)
        self.assertEqual(DEFAULT_TRIG_CLASSES, params.trigPointClasses)

    def testDistanceBaseZeroValueError(self):
        """
        If the distanceBase is zero, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        self.assertRaisesRegex(ValueError, error, Parameters, [], [],
                               distanceBase=0.0)

    def testDistanceBaseLessThanZeroValueError(self):
        """
        If the distanceBase is < 0, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        self.assertRaisesRegex(ValueError, error, Parameters, [], [],
                               distanceBase=-1.0)

    def testSaveReturnsItsArgument(self):
        """
        The save function must return its (fp) argument.
        """
        params = Parameters([AlphaHelix], [Peaks], limitPerLandmark=3,
                            maxDistance=10)
        fp = StringIO()
        self.assertIs(fp, params.save(fp))

    def testSave(self):
        """
        Saving parameters as JSON must work correctly.
        """
        params = Parameters([AlphaHelix], [Peaks], limitPerLandmark=3,
                            maxDistance=19, minDistance=5, distanceBase=1.2)

        fp = StringIO()
        params.save(fp)
        expected = {
            'landmarkClassNames': ['AlphaHelix'],
            'trigPointClassNames': ['Peaks'],
            'limitPerLandmark': 3,
            'maxDistance': 19,
            'minDistance': 5,
            'distanceBase': 1.2,
        }
        self.assertEqual(expected, loads(fp.getvalue()))

    def testRestore(self):
        """
        The restore method must produce the same parameter values that were
        present in a Parameters instance when it is saved.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            distanceBase=3.0, limitPerLandmark=10,
                            maxDistance=77, minDistance=66)
        fp = StringIO()
        params.save(fp)
        fp.seek(0)
        result = Parameters.restore(fp)

        self.assertEqual([AlphaHelix, BetaStrand], result.landmarkClasses)
        self.assertEqual([Peaks, Troughs], result.trigPointClasses)
        self.assertEqual(3.0, result.distanceBase)
        self.assertEqual(10, result.limitPerLandmark)
        self.assertEqual(77, result.maxDistance)
        self.assertEqual(66, result.minDistance)

    def testFinderOrderInvariant(self):
        """
        The parameter checksum must be identical when finders are specified,
        no matter what order the finders are given.
        """
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        params2 = Parameters([BetaStrand, AlphaHelix], [Troughs, Peaks])
        self.assertEqual(params1.checksum, params2.checksum)

    def testSaveLoadMissingLandmark(self):
        """
        If parameters are saved with a landmark class whose name cannot be
        found when the parameters are later loaded, a ValueError error must be
        raised.
        """
        params = Parameters([AlphaHelix], [])
        fp = StringIO()
        params.save(fp)
        newSave = fp.getvalue().replace('AlphaHelix', 'Non-existent')
        error = ('^Could not find landscape finder class Non-existent. '
                 'Has that class been renamed or removed\?$')
        self.assertRaisesRegex(ValueError, error, Parameters.restore,
                               StringIO(newSave))

    def testSaveLoadMissingTrigPoint(self):
        """
        If parameters are saved with a trig point class whose name cannot be
        found when the parameters are later loaded, a ValueError error must be
        raised.
        """
        params = Parameters([], [Peaks])
        fp = StringIO()
        params.save(fp)
        newSave = fp.getvalue().replace('Peaks', 'Non-existent')
        error = ('^Could not find trig point finder class Non-existent. '
                 'Has that class been renamed or removed\?$')
        self.assertRaisesRegex(ValueError, error, Parameters.restore,
                               StringIO(newSave))

    def testRestoreInvalidJSON(self):
        """
        If parameter restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        error = '^Expected object or value$'
        self.assertRaisesRegex(ValueError, error, Parameters.restore,
                               StringIO('xxx'))

    def testCompareIdentical(self):
        """
        The compare method must return None when two parameter instances
        have the same values.
        """
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                             distanceBase=3.0, limitPerLandmark=10,
                             maxDistance=77, minDistance=66)
        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                             distanceBase=3.0, limitPerLandmark=10,
                             maxDistance=77, minDistance=66)
        self.assertIs(None, params1.compare(params2))

    def testCompareDifferentLandmarkFinders(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different landmark finders.
        """
        params1 = Parameters([AlphaHelix, BetaStrand], [])
        params2 = Parameters([AlphaHelix], [])
        expected = ("Summary of differences:\n"
                    "\tParam 'landmarkClasses' values "
                    "[<class 'light.landmarks.alpha_helix.AlphaHelix'>, "
                    "<class 'light.landmarks.beta_strand.BetaStrand'>] and "
                    "[<class 'light.landmarks.alpha_helix.AlphaHelix'>] "
                    "differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentTrigPointFinders(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different trig point finders.
        """
        params1 = Parameters([], [Peaks, Troughs])
        params2 = Parameters([], [Peaks])
        expected = ("Summary of differences:\n"
                    "\tParam 'trigPointClasses' values "
                    "[<class 'light.trig.peaks.Peaks'>, "
                    "<class 'light.trig.troughs.Troughs'>] and "
                    "[<class 'light.trig.peaks.Peaks'>] "
                    "differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentDistanceBase(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different distance bases.
        """
        params1 = Parameters([], [], distanceBase=1.0)
        params2 = Parameters([], [], distanceBase=2.0)
        expected = ("Summary of differences:\n"
                    "\tParam 'distanceBase' values 1.0 and 2.0 differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentLimitPerLandmark(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different limit per landmark.
        """
        params1 = Parameters([], [], limitPerLandmark=10)
        params2 = Parameters([], [], limitPerLandmark=20)
        expected = ("Summary of differences:\n"
                    "\tParam 'limitPerLandmark' values 10 and 20 differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentMaxDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different max distances.
        """
        params1 = Parameters([], [], maxDistance=10)
        params2 = Parameters([], [], maxDistance=20)
        expected = ("Summary of differences:\n"
                    "\tParam 'maxDistance' values 10 and 20 differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentMinDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different min distances.
        """
        params1 = Parameters([], [], minDistance=10)
        params2 = Parameters([], [], minDistance=20)
        expected = ("Summary of differences:\n"
                    "\tParam 'minDistance' values 10 and 20 differ.")
        self.assertEqual(expected, params1.compare(params2))

    def testCompareDifferentMaxDistanceAndMinDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have multiple different attributes.
        """
        params1 = Parameters([], [], maxDistance=10, minDistance=30)
        params2 = Parameters([], [], maxDistance=20, minDistance=40)
        expected = ("Summary of differences:\n"
                    "\tParam 'maxDistance' values 10 and 20 differ.\n"
                    "\tParam 'minDistance' values 30 and 40 differ.")
        self.assertEqual(expected, params1.compare(params2))
