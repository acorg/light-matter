from unittest import TestCase
from json import loads
from io import StringIO

from light.parameters import Parameters
from light.landmarks import DEFAULT_LANDMARK_CLASSES, AlphaHelix, BetaStrand
from light.trig import DEFAULT_TRIG_CLASSES, Peaks, Troughs


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
