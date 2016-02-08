import six
import argparse
from unittest import TestCase
from json import loads
from six import StringIO

from light.parameters import DatabaseParameters, FindParameters
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
            '--deltaScale', '0.2',
        ])

        # Parsing must do the expected thing.
        self.assertEqual('Always', args.significanceMethod)
        self.assertEqual(0.4, args.significanceFraction)
        self.assertEqual('MinHashesScore', args.scoreMethod)
        self.assertEqual(5, args.featureMatchScore)
        self.assertEqual(6, args.featureMismatchScore)
        self.assertEqual(0.2, args.deltaScale)

        # We must be able to make an instance from the parsed args.
        findParams = FindParameters.fromArgs(args)
        self.assertEqual('Always', findParams.significanceMethod)
        self.assertEqual(0.4, findParams.significanceFraction)
        self.assertEqual('MinHashesScore', findParams.scoreMethod)
        self.assertEqual(5, findParams.featureMatchScore)
        self.assertEqual(6, findParams.featureMismatchScore)
        self.assertEqual(0.2, findParams.deltaScale)


class TestDatabaseParameters(TestCase):
    """
    Tests for the light.database.DatabaseParameters class.
    """
    def testDefaults(self):
        """
        If no specific parameter values are given, the defaults must be set.
        """
        dbParams = DatabaseParameters()
        self.assertEqual([cls.NAME for cls in DEFAULT_LANDMARK_CLASSES],
                         dbParams.landmarkFinderNames())
        self.assertEqual([cls.NAME for cls in DEFAULT_TRIG_CLASSES],
                         dbParams.trigPointFinderNames())
        self.assertEqual(DatabaseParameters.DEFAULT_DISTANCE_BASE,
                         dbParams.distanceBase)
        self.assertEqual(DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE,
                         dbParams.featureLengthBase)
        self.assertEqual(DatabaseParameters.DEFAULT_LIMIT_PER_LANDMARK,
                         dbParams.limitPerLandmark)
        self.assertEqual(DatabaseParameters.DEFAULT_MAX_DISTANCE,
                         dbParams.maxDistance)
        self.assertEqual(DatabaseParameters.DEFAULT_MIN_DISTANCE,
                         dbParams.minDistance)
        self.assertEqual(DatabaseParameters.DEFAULT_RANDOM_LANDMARK_DENSITY,
                         dbParams.randomLandmarkDensity)
        self.assertEqual(DatabaseParameters.DEFAULT_RANDOM_TRIG_POINT_DENSITY,
                         dbParams.randomTrigPointDensity)

    def testNotDefaults(self):
        """
        If specific parameter values are given, the passed values must be set.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                      trigPoints=[Peaks, Troughs],
                                      distanceBase=3.0, featureLengthBase=1.7,
                                      limitPerLandmark=10, maxDistance=77,
                                      minDistance=66,
                                      randomLandmarkDensity=0.1,
                                      randomTrigPointDensity=0.2)
        self.assertEqual(['AlphaHelix', 'BetaStrand'],
                         dbParams.landmarkFinderNames())
        self.assertEqual(['Peaks', 'Troughs'], dbParams.trigPointFinderNames())
        self.assertEqual(3.0, dbParams.distanceBase)
        self.assertEqual(1.7, dbParams.featureLengthBase)
        self.assertEqual(10, dbParams.limitPerLandmark)
        self.assertEqual(77, dbParams.maxDistance)
        self.assertEqual(66, dbParams.minDistance)
        self.assertEqual(0.1, dbParams.randomLandmarkDensity)
        self.assertEqual(0.2, dbParams.randomTrigPointDensity)

    def testDistanceBaseZeroValueError(self):
        """
        If the distanceBase is zero, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        six.assertRaisesRegex(self, ValueError, error, DatabaseParameters,
                              distanceBase=0.0)

    def testDistanceBaseLessThanZeroValueError(self):
        """
        If the distanceBase is < 0, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        six.assertRaisesRegex(self, ValueError, error, DatabaseParameters,
                              distanceBase=-1.0)

    def testFeatureLengthBaseZeroValueError(self):
        """
        If the featureLengthBase is zero, a ValueError must be raised.
        """
        error = 'featureLengthBase must be > 0\\.'
        six.assertRaisesRegex(self, ValueError, error, DatabaseParameters,
                              featureLengthBase=0.0)

    def testFeatureLengthBaseLessThanZeroValueError(self):
        """
        If the featureLengthBase is < 0, a ValueError must be raised.
        """
        error = 'featureLengthBase must be > 0\\.'
        six.assertRaisesRegex(self, ValueError, error, DatabaseParameters,
                              featureLengthBase=-1.0)

    def testSaveReturnsItsArgument(self):
        """
        The save function must return its (fp) argument.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix],
                                      trigPoints=[Peaks], limitPerLandmark=3,
                                      maxDistance=10)
        fp = StringIO()
        self.assertIs(fp, dbParams.save(fp))

    def testSave(self):
        """
        Saving parameters as JSON must work correctly.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix],
                                      trigPoints=[Peaks], limitPerLandmark=3,
                                      maxDistance=19, minDistance=5,
                                      distanceBase=1.2, featureLengthBase=1.7)
        fp = StringIO()
        dbParams.save(fp)
        expected = {
            'landmarks': ['AlphaHelix'],
            'trigPoints': ['Peaks'],
            'limitPerLandmark': 3,
            'maxDistance': 19,
            'minDistance': 5,
            'distanceBase': 1.2,
            'featureLengthBase': 1.7,
            'randomLandmarkDensity': 0.1,
            'randomTrigPointDensity': 0.1,
        }
        self.assertEqual(expected, loads(fp.getvalue()))

    def testRestore(self):
        """
        The restore method must produce the same parameter values that were
        present in a Parameters instance when it is saved.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                      trigPoints=[Peaks, Troughs],
                                      distanceBase=3.0, featureLengthBase=3.9,
                                      limitPerLandmark=10, maxDistance=77,
                                      minDistance=66,
                                      randomLandmarkDensity=0.1,
                                      randomTrigPointDensity=0.9)
        fp = StringIO()
        dbParams.save(fp)
        fp.seek(0)
        result = DatabaseParameters.restore(fp)

        self.assertEqual([AlphaHelix(), BetaStrand()], result.landmarkFinders)
        self.assertEqual([Peaks(), Troughs()], result.trigPointFinders)
        self.assertEqual(3.0, result.distanceBase)
        self.assertEqual(3.9, result.featureLengthBase)
        self.assertEqual(10, result.limitPerLandmark)
        self.assertEqual(77, result.maxDistance)
        self.assertEqual(66, result.minDistance)
        self.assertEqual(0.1, result.randomLandmarkDensity)
        self.assertEqual(0.9, result.randomTrigPointDensity)

    def testFinderOrderInvariant(self):
        """
        The parameter checksum must be identical when finders are specified,
        no matter what order the finders are given.
        """
        dbParams1 = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                       trigPoints=[Peaks, Troughs])
        dbParams2 = DatabaseParameters(landmarks=[BetaStrand, AlphaHelix],
                                       trigPoints=[Troughs, Peaks])
        self.assertEqual(dbParams1.checksum, dbParams2.checksum)

    def testSaveLoadMissingLandmark(self):
        """
        If parameters are saved with a landmark class whose name cannot be
        found when the parameters are later loaded, a ValueError error must be
        raised.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        fp = StringIO()
        dbParams.save(fp)
        newSave = fp.getvalue().replace('AlphaHelix', 'Non-existent')
        error = '^Could not find landmark finder class \'Non-existent\'\.$'
        six.assertRaisesRegex(self, ValueError, error,
                              DatabaseParameters.restore, StringIO(newSave))

    def testSaveLoadMissingTrigPoint(self):
        """
        If parameters are saved with a trig point class whose name cannot be
        found when the parameters are later loaded, a ValueError error must be
        raised.
        """
        dbParams = DatabaseParameters(trigPoints=[Peaks])
        fp = StringIO()
        dbParams.save(fp)
        newSave = fp.getvalue().replace('Peaks', 'Non-existent')
        error = '^Could not find trig point finder class \'Non-existent\'.$'
        six.assertRaisesRegex(self, ValueError, error,
                              DatabaseParameters.restore, StringIO(newSave))

    def testRestoreInvalidJSON(self):
        """
        If parameter restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        error = '^Expected object or value$'
        six.assertRaisesRegex(self, ValueError, error,
                              DatabaseParameters.restore, StringIO('xxx'))

    def testCompareIdentical(self):
        """
        The compare method must return None when two parameter instances
        have the same values.
        """
        dbParams1 = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                       trigPoints=[Peaks, Troughs],
                                       distanceBase=3.0, featureLengthBase=3.9,
                                       limitPerLandmark=10, maxDistance=77,
                                       minDistance=66,
                                       randomLandmarkDensity=0.1,
                                       randomTrigPointDensity=0.2)
        dbParams2 = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                       trigPoints=[Peaks, Troughs],
                                       distanceBase=3.0, featureLengthBase=3.9,
                                       limitPerLandmark=10, maxDistance=77,
                                       minDistance=66,
                                       randomLandmarkDensity=0.1,
                                       randomTrigPointDensity=0.2)
        self.assertIs(None, dbParams1.compare(dbParams2))

    def testCompareDifferentLandmarkFinders(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different landmark finders.
        """
        dbParams1 = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand])
        dbParams2 = DatabaseParameters(landmarks=[AlphaHelix])
        expected = ("Summary of differences:\n"
                    "  Param 'landmarks' values "
                    "['AlphaHelix', 'BetaStrand'] and ['AlphaHelix'] differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentTrigPointFinders(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different trig point finders.
        """
        dbParams1 = DatabaseParameters(trigPoints=[Peaks, Troughs])
        dbParams2 = DatabaseParameters(trigPoints=[Peaks])
        expected = ("Summary of differences:\n"
                    "  Param 'trigPoints' values "
                    "['Peaks', 'Troughs'] and ['Peaks'] differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentDistanceBase(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different distance bases.
        """
        dbParams1 = DatabaseParameters(distanceBase=1.0)
        dbParams2 = DatabaseParameters(distanceBase=2.0)
        expected = ("Summary of differences:\n"
                    "  Param 'distanceBase' values 1.0 and 2.0 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentFeatureLengthBase(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different featureLength bases.
        """
        dbParams1 = DatabaseParameters(featureLengthBase=1.0)
        dbParams2 = DatabaseParameters(featureLengthBase=2.0)
        expected = ("Summary of differences:\n"
                    "  Param 'featureLengthBase' values 1.0 and 2.0 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentLimitPerLandmark(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different limit per landmark.
        """
        dbParams1 = DatabaseParameters(limitPerLandmark=10)
        dbParams2 = DatabaseParameters(limitPerLandmark=20)
        expected = ("Summary of differences:\n"
                    "  Param 'limitPerLandmark' values 10 and 20 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentMaxDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different max distances.
        """
        dbParams1 = DatabaseParameters(maxDistance=10)
        dbParams2 = DatabaseParameters(maxDistance=20)
        expected = ("Summary of differences:\n"
                    "  Param 'maxDistance' values 10 and 20 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentMinDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different min distances.
        """
        dbParams1 = DatabaseParameters(minDistance=10)
        dbParams2 = DatabaseParameters(minDistance=20)
        expected = ("Summary of differences:\n"
                    "  Param 'minDistance' values 10 and 20 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentMaxDistanceAndMinDistance(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have multiple different attributes.
        """
        dbParams1 = DatabaseParameters(maxDistance=10, minDistance=30)
        dbParams2 = DatabaseParameters(maxDistance=20, minDistance=40)
        expected = ("Summary of differences:\n"
                    "  Param 'maxDistance' values 10 and 20 differ.\n"
                    "  Param 'minDistance' values 30 and 40 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testLandmarkFinderNamesNoLandmarkFinders(self):
        """
        The landmarkFinderNames method must return an empty list when there
        are no landmark finders.
        """
        dbParams = DatabaseParameters(landmarks=[])
        self.assertEqual([], dbParams.landmarkFinderNames())

    def testLandmarkFinderNames(self):
        """
        The landmarkFinderNames method must return the expected list of names
        of landmark finders.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand])
        self.assertEqual(['AlphaHelix', 'BetaStrand'],
                         dbParams.landmarkFinderNames())

    def testTrigPointFinderNamesNoTrigPointFinders(self):
        """
        The trigPointFinderNames method must return an empty list when there
        are no trig point finders.
        """
        dbParams = DatabaseParameters(trigPoints=[])
        self.assertEqual([], dbParams.trigPointFinderNames())

    def testTrigPointFinderNames(self):
        """
        The trigPointFinderNames method must return the expected list of names
        of trig point finders.
        """
        dbParams = DatabaseParameters(trigPoints=[Peaks, Troughs])
        self.assertEqual(['Peaks', 'Troughs'],
                         dbParams.trigPointFinderNames())
