import six
import argparse
from unittest import TestCase, skip
from json import loads
from six import StringIO

from light.bin_score import MinHashesScore
from light.landmarks import DEFAULT_LANDMARK_CLASSES, AlphaHelix, BetaStrand
from light.overall_score import BestBinScore
from light.parameters import DatabaseParameters, FindParameters
from light.significance import HashFraction
from light.trig import DEFAULT_TRIG_CLASSES, Peaks, Troughs


class TestFindParameters(TestCase):
    """
    Tests for the light.database.FindParameters class.
    """
    def testDefaults(self):
        """
        If no specific parameter values are given, the defaults must be set.
        """
        findParams = FindParameters()
        self.assertEqual(FindParameters.DEFAULT_SIGNIFICANCE_METHOD,
                         findParams.significanceMethod)
        self.assertEqual(FindParameters.DEFAULT_SCORE_METHOD,
                         findParams.scoreMethod)
        self.assertEqual(FindParameters.DEFAULT_OVERALL_SCORE_METHOD,
                         findParams.overallScoreMethod)
        self.assertEqual(FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         findParams.featureMatchScore)
        self.assertEqual(FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
                         findParams.featureMismatchScore)
        self.assertEqual(FindParameters.DEFAULT_WEIGHTS,
                         findParams.weights)
        self.assertEqual(FindParameters.DEFAULT_DELTA_SCALE,
                         findParams.deltaScale)

    def testDefaultValues(self):
        """
        Guard against accidental edits of default parameter values.
        """
        self.assertEqual(HashFraction.__name__,
                         FindParameters.DEFAULT_SIGNIFICANCE_METHOD)
        self.assertEqual(MinHashesScore.__name__,
                         FindParameters.DEFAULT_SCORE_METHOD)
        self.assertEqual(BestBinScore.__name__,
                         FindParameters.DEFAULT_OVERALL_SCORE_METHOD)
        self.assertEqual(1.0, FindParameters.DEFAULT_FEATURE_MATCH_SCORE)
        self.assertEqual(-1.0, FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE)
        self.assertEqual(
            {
                'AlphaHelix': 1.0,
                'AlphaHelix_3_10': 1.0,
                'AlphaHelix_pi': 1.0,
                'BetaStrand': 1.0,
                'BetaTurn': 1.0,
                'AminoAcidsLm': 1.0,
                'GOR4AlphaHelix': 1.0,
                'GOR4BetaStrand': 1.0,
                'GOR4Coil': 1.0,
                'Prosite': 1.0,
                'Peaks': 1.0,
                'Troughs': 1.0,
                'AminoAcids': 1.0,
                'IndividualPeaks': 1.0,
                'IndividualTroughs': 1.0,
            },
            FindParameters.DEFAULT_WEIGHTS)
        self.assertEqual(1.0, FindParameters.DEFAULT_DELTA_SCALE)

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

    def testUnknownMethods(self):
        """
        It is currently not possible to test passing unknown methods.
        """
        # parse_args prints to stderr and calls sys.exit if an argument
        # whose possible choices are not met is encountered. It also
        # imports sys before we get a chance to patch sys.exit.  I'm
        # leaving this (non-)test here so you can see I (Terry) tried to
        # treat these cases and also because it may become useful if
        # argparse becomes more flexible.
        #
        # This problem exists for testing the 'choices' arguments:
        # significanceMethod, scoreMethod, overallScoreMethod.


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

    def testDefaultValues(self):
        """
        Guard against accidental edits of default parameter values.
        """
        # DEFAULT_LANDMARK_CLASSES and DEFAULT_TRIG_CLASSES are tested in
        # test/landmarks/test_landmark.py and test/trig/test_trig.py
        self.assertEqual(10, DatabaseParameters.DEFAULT_LIMIT_PER_LANDMARK)
        self.assertEqual(200, DatabaseParameters.DEFAULT_MAX_DISTANCE)
        self.assertEqual(1, DatabaseParameters.DEFAULT_MIN_DISTANCE)
        self.assertEqual(1.1, DatabaseParameters.DEFAULT_DISTANCE_BASE)
        self.assertEqual(1.35, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        self.assertEqual(0.1,
                         DatabaseParameters.DEFAULT_RANDOM_LANDMARK_DENSITY)
        self.assertEqual(0.1,
                         DatabaseParameters.DEFAULT_RANDOM_TRIG_POINT_DENSITY)

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
        dbParams = DatabaseParameters(landmarks=[], trigPoints=[])
        fp = StringIO()
        self.assertIs(fp, dbParams.save(fp))

    def testSave(self):
        """
        Saving parameters as JSON must work correctly.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix],
                                      trigPoints=[Peaks], limitPerLandmark=3,
                                      maxDistance=19, minDistance=5,
                                      distanceBase=1.2, featureLengthBase=1.7,
                                      randomLandmarkDensity=0.3,
                                      randomTrigPointDensity=0.9)
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
            'randomLandmarkDensity': 0.3,
            'randomTrigPointDensity': 0.9,
        }
        self.assertEqual(expected, loads(fp.getvalue()))

    def testRestore(self):
        """
        The restore method must produce the same parameter values that were
        present in a DatabaseParameters instance when it is saved.
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

    def testCompareDifferentRandomLandmarkDensity(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different randomLandmarkDensity
        values.
        """
        dbParams1 = DatabaseParameters(randomLandmarkDensity=0.2)
        dbParams2 = DatabaseParameters(randomLandmarkDensity=0.4)
        expected = (
            "Summary of differences:\n"
            "  Param 'randomLandmarkDensity' values 0.2 and 0.4 differ.")
        self.assertEqual(expected, dbParams1.compare(dbParams2))

    def testCompareDifferentRandomTrigPointDensity(self):
        """
        The compare method must return a description of parameter differences
        when two parameter instances have different randomTrigPointDensity
        values.
        """
        dbParams1 = DatabaseParameters(randomTrigPointDensity=0.2)
        dbParams2 = DatabaseParameters(randomTrigPointDensity=0.4)
        expected = (
            "Summary of differences:\n"
            "  Param 'randomTrigPointDensity' values 0.2 and 0.4 differ.")
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

    def testFromArgsNoArgs(self):
        """
        If no arguments are given on a command line, default parameters must
        be returned by fromArgs.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args([])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertIs(None, dbParams.compare(DatabaseParameters()))

    def testFromArgsNoLandmarks(self):
        """
        If --noLandmarks is given on a command line, the returned parameters
        from fromArgs must have no landmarks.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args(['--noLandmarks'])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertEqual([], dbParams.landmarkFinders)

    @skip('Cannot test errors causing argparse to call sys.exit')
    def testFromArgsNonexistentLandmark(self):
        """
        If --landmark is given on a command line with a non-existent landmark
        name, we should be able to catch it but currently cannot.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        # The following doesn't work, as parse_args prints to stderr and
        # calls sys.exit if an argument whose possible choices are not met
        # is encountered. It also imports sys before we get a chance to
        # patch sys.exit.  I'm leaving this (non-)test here so you can see
        # I (Terry) tried to treat this case and also because it may become
        # useful if argparse becomes more flexible.
        #
        # args = parser.parse_args(['--landmark', 'non-existent'])

    def testFromArgsLandmarks(self):
        """
        If --landmark is given on a command line, the returned parameters
        from fromArgs must have the expected landmark finders.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args(
            ['--landmark', 'AlphaHelix', '--landmark', 'BetaStrand'])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertEqual(['AlphaHelix', 'BetaStrand'],
                         dbParams.landmarkFinderNames())

    def testFromArgsNoTrigPoints(self):
        """
        If --noTrigPoints is given on a command line, the returned parameters
        from fromArgs must have no trig points.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args(['--noTrigPoints'])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertEqual([], dbParams.trigPointFinders)

    @skip('Cannot test errors causing argparse to call sys.exit')
    def testFromArgsNonexistentTrigPoint(self):
        """
        If --trig is given on a command line with a non-existent trig point
        name, we should be able to catch it but currently cannot.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        # The following doesn't work, as parse_args prints to stderr and
        # calls sys.exit if an argument whose possible choices are not met
        # is encountered. It also imports sys before we get a chance to
        # patch sys.exit.  I'm leaving this (non-)test here so you can see
        # I (Terry) tried to treat this case and also because it may become
        # useful if argparse becomes more flexible.
        #
        # args = parser.parse_args(['--trig', 'non-existent'])

    def testFromArgsTrigPoints(self):
        """
        If --trig is given on a command line, the returned parameters
        from fromArgs must have the expected trig point finders.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args(['--trig', 'Troughs', '--trig', 'Peaks'])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertEqual(['Peaks', 'Troughs'],
                         sorted(dbParams.trigPointFinderNames()))

    def testFromArgsScalarParameters(self):
        """
        All scalar arguments must be processed by fromArgs as expected.
        """
        parser = argparse.ArgumentParser()
        DatabaseParameters.addArgsToParser(parser)
        args = parser.parse_args([
            '--limitPerLandmark', '5',
            '--maxDistance', '10',
            '--minDistance', '3',
            '--distanceBase', '1.9',
            '--featureLengthBase', '2.3',
            '--randomLandmarkDensity', '0.7',
            '--randomTrigPointDensity', '0.3',
        ])
        dbParams = DatabaseParameters.fromArgs(args)
        self.assertEqual(5, dbParams.limitPerLandmark)
        self.assertEqual(10, dbParams.maxDistance)
        self.assertEqual(3, dbParams.minDistance)
        self.assertEqual(1.9, dbParams.distanceBase)
        self.assertEqual(2.3, dbParams.featureLengthBase)
        self.assertEqual(0.7, dbParams.randomLandmarkDensity)
        self.assertEqual(0.3, dbParams.randomTrigPointDensity)
