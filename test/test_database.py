from unittest import TestCase
from json import loads
from cStringIO import StringIO
from hashlib import sha256

from light.features import Landmark, TrigPoint
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.beta_strand import BetaStrand
from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.database import Database
from dark.reads import AARead


class TestDatabase(TestCase):
    """
    Tests for the light.database.Database class.
    """
    @staticmethod
    def _update(s, checksum):
        """
        Add the string representation of an object to a checksum,
        followed by a NUL.

        @param s: Anything that can be converted to a C{str} to be added to
            the checksum.
        @param checksum: Any hashlib sum instance with an update method.
        """
        checksum.update(str(s) + '\0')

    def testAboveMeanThresholdDefault(self):
        """
        The above mean threshold default value must be as expected.
        """
        self.assertEqual(15, Database.ABOVE_MEAN_THRESHOLD_DEFAULT)

    def testFindersAreStored(self):
        """
        The list of landmark and trig point finders must be stored correctly.
        """
        db = Database([AlphaHelix], [Peaks])
        self.assertEqual([AlphaHelix], db.landmarkFinderClasses)
        self.assertEqual([Peaks], db.trigPointFinderClasses)

    def testInitialStatistics(self):
        """
        The database statistics must be initially correct.
        """
        db = Database([], [])
        self.assertEqual(0, db.subjectCount)
        self.assertEqual(0, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testInitialDatabaseIsEmpty(self):
        """
        The database must be empty if no reads have been added.
        """
        db = Database([AlphaHelix], [Peaks])
        self.assertEqual({}, db.d)

    def testKey(self):
        """
        The database key function must return the expected value.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        self.assertEqual('A1:B:10', db.key(landmark, trigPoint))

    def testKeyWithRepeatCount(self):
        """
        The database key function must return the expected value when the
        landmark it is passed has a repeat count.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0, 5)
        trigPoint = TrigPoint('name', 'B', 10)
        self.assertEqual('A5:B:10', db.key(landmark, trigPoint))

    def testInitialDatabaseHasNoReadInfo(self):
        """
        The database must not have any stored read information if no reads have
        been added.
        """
        db = Database([], [])
        self.assertEqual([], db.subjectInfo)

    def testOneReadOneLandmark(self):
        """
        If one read is added but it only has one landmark, nothing is added
        to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkReadInfo(self):
        """
        If one read is added an entry is appended to the read info.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual([('id', 'FRRRFRRRF')], db.subjectInfo)

    def testOneReadOneLandmarkStatistics(self):
        """
        If one read is added the database statistics must be correct.
        """
        db = Database([], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(9, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testOneReadTwoLandmarks(self):
        """
        If one read is added and it has two landmarks, two keys are added
        to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual({'A3:A2:23': [{'subjectIndex': 0,
                                        'offset': 23}],
                         'A2:A3:-23': [{'subjectIndex': 0,
                                        'offset': 0}]},
                         db.d)

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one read is added, the database statistics must be correct.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(36, db.totalResidues)
        self.assertEqual(22, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarks(self):
        """
        If two identical reads are added, both with two landmarks, two keys
        are added to the database and both reads are listed in the dictionary
        values for those two keys.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual({'A3:A2:23': [{'subjectIndex': 0,
                                        'offset': 23},
                                       {'subjectIndex': 1,
                                        'offset': 23}],
                         'A2:A3:-23': [{'subjectIndex': 0,
                                        'offset': 0},
                                       {'subjectIndex': 1,
                                        'offset': 0}]},
                         db.d)

    def testTwoReadsTwoLandmarksStatistics(self):
        """
        If two identical reads are added, the database statistics must be
        correct.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual(2, db.subjectCount)
        self.assertEqual(72, db.totalResidues)
        self.assertEqual(44, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarksLimitZeroPairsPerLandmark(self):
        """
        If two identical reads are added, both with two landmarks, no keys
        will be added to the dictionary if limitPerLandmark is zero.
        """
        db = Database([AlphaHelix], [], limitPerLandmark=0)
        db.addSubject(AARead('id1', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual({}, db.d)

    def testTwoReadsTwoLandmarksDifferentOffsets(self):
        """
        If two reads are added, both with two landmarks separated by the same
        distance, two keys are added to the database and both reads are listed
        in the dictionary values for those two keys.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'AFRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        self.assertEqual({'A3:A2:23': [{'subjectIndex': 0,
                                        'offset': 24},
                                       {'subjectIndex': 1,
                                        'offset': 23}],
                          'A2:A3:-23': [{'subjectIndex': 0,
                                         'offset': 1},
                                        {'subjectIndex': 1,
                                         'offset': 0}]},
                         db.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one read is added and it has one landmark and one peak, one pair is
        added to the database.
        """
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual({'A2:P:-10': [{'subjectIndex': 0,
                                        'offset': 0}]},
                         db.d)

    def testOneReadOneLandmarkOnePeakNoTrigFinders(self):
        """
        If one read is added and it has one landmark and one peak, but no trig
        finders are in use, nothing is added to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkTwoPeaks(self):
        """
        If one read is added and it has one landmark and two peaks, two pairs
        are added to the database.
        """
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({'A2:P:-13': [{'subjectIndex': 0,
                                        'offset': 0}],
                          'A2:P:-10': [{'subjectIndex': 0,
                                        'offset': 0}]},
                         db.d)

    def testOneReadOneLandmarkTwoPeaksLimitOnePairPerLandmark(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        limit of one pair per landmarks is imposed, only one key is added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], limitPerLandmark=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({'A2:P:-10': [{'subjectIndex': 0,
                                        'offset': 0}]},
                         db.d)

    def testOneReadOneLandmarkTwoPeaksSevereMaxDistance(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        severe maximum distance is imposed, no keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMaxDistance(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        maximum distance is imposed that makes one of the peaks too far
        away, only one key is added to the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({'A2:P:-10': [{'subjectIndex': 0,
                                        'offset': 0}]},
                         db.d)

    def testOneReadOneLandmarkTwoPeaksLargeMaxDistance(self):
        """
        If one read is added and it has one landmark and two peaks, and a
        maximum distance is imposed that is greater than the distance to the
        peaks, two keys are added to the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=15)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({'A2:P:-13': [{'subjectIndex': 0,
                                        'offset': 0}],
                         'A2:P:-10': [{'subjectIndex': 0,
                                       'offset': 0}]},
                         db.d)

    def testSaveLoadEmpty(self):
        """
        When asked to save and then load an empty database, the correct
        database must result.
        """
        db = Database([], [])
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.load(fp)
        self.assertEqual(0, result.subjectCount)
        self.assertEqual(0, result.totalCoveredResidues)
        self.assertEqual({}, result.d)
        self.assertEqual([], result.landmarkFinderClasses)
        self.assertEqual(None, result.limitPerLandmark)
        self.assertEqual(None, result.maxDistance)
        self.assertEqual([], result.trigPointFinderClasses)
        self.assertEqual(0, result.totalResidues)

    def testSaveLoadNonEmpty(self):
        """
        When asked to save and then load a non-empty database, the correct
        database must result.
        """
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.load(fp)
        self.assertEqual(db.subjectCount, result.subjectCount)
        self.assertEqual(db.totalCoveredResidues, result.totalCoveredResidues)
        self.assertEqual(db.d, result.d)
        self.assertEqual(db.landmarkFinderClasses,
                         result.landmarkFinderClasses)
        self.assertEqual(db.limitPerLandmark, result.limitPerLandmark)
        self.assertEqual(db.maxDistance, result.maxDistance)
        self.assertEqual(db.trigPointFinderClasses,
                         result.trigPointFinderClasses)
        self.assertEqual(db.totalResidues, result.totalResidues)

    def testFindNotMatching(self):
        """
        A non-matching key must not be found.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRR')
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindOneMatchingInsignificant(self):
        """
        One matching subject should be found, but is not significant with the
        default value of aboveMeanThreshold.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual(
            {
                0: [
                    {
                        'trigPointName': 'Peaks',
                        'distance': 10,
                        'landmarkLength': 9,
                        'readOffset': 0,
                        'subjectOffset': 1,
                        'landmarkName': 'AlphaHelix',
                    },
                ],
            },
            result.matches)
        self.assertEqual(0, len(list(result.significant())))

    def testFindOneMatchingSignificant(self):
        """
        One matching and significant subject must be found if the
        aboveMeanThreshold is sufficiently low.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(subject)
        result = db.find(query, aboveMeanThreshold=0.1)
        self.assertEqual(
            {
                0: [
                    {
                        'trigPointName': 'Peaks',
                        'distance': 10,
                        'landmarkLength': 9,
                        'readOffset': 0,
                        'subjectOffset': 1,
                        'landmarkName': 'AlphaHelix',
                    },
                ],
            },
            result.matches)

    def testFindNoneMatchingTooSmallDistance(self):
        """
        One matching key must be found.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks], maxDistance=1)
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindNoneMatchingNoTrigPoint(self):
        """
        One matching key must be found.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindTwoMatchingInSameSubject(self):
        """
        Two matching keys in the same subject must be found.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual(
            {
                0: [
                    {
                        'distance': 10,
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'readOffset': 0,
                        'subjectOffset': 0,
                        'trigPointName': 'Peaks',
                    },
                    {
                        'distance': 13,
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'readOffset': 0,
                        'subjectOffset': 0,
                        'trigPointName': 'Peaks',
                    }
                ],
            },
            result.matches)

    def testSaveParamsAsJSONReturnsItsArgument(self):
        """
        The saveParamsAsJSON function must return its (fp) argument.
        """
        db = Database([AlphaHelix], [Peaks], limitPerLandmark=3,
                      maxDistance=10)
        io = StringIO()
        self.assertIs(io, db.saveParamsAsJSON(io))

    def testSaveParamsAsJSON(self):
        """
        Saving parameters as JSON must work correctly.
        """
        db = Database([AlphaHelix], [Peaks], limitPerLandmark=3,
                      maxDistance=10)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))

        checksum = sha256()
        self._update(AlphaHelix.NAME, checksum)
        self._update(Peaks.NAME, checksum)
        self._update(3, checksum)  # Limit per landmark value.
        self._update(10, checksum)  # Max distance value.
        self._update('id', checksum)
        self._update('FRRRFRRRFASAASA', checksum)

        # The database dictionary now contains
        #
        # {
        #     'A2:P:-10': [{'subjectIndex': 0, 'offset': 0}]}
        # }
        #
        # We have to add the keys in order and the keys and values of the
        # values in order.

        self._update('A2:P:-10', checksum)
        self._update('offset', checksum)
        self._update(0, checksum)
        self._update('subjectIndex', checksum)
        self._update(0, checksum)

        out = StringIO()
        db.saveParamsAsJSON(out)
        expected = {
            'checksum': checksum.hexdigest(),
            'landmarkFinderClasses': ['AlphaHelix'],
            'trigPointFinderClasses': ['Peaks'],
            'limitPerLandmark': 3,
            'maxDistance': 10,
            'subjectCount': 1,
            'totalResidues': 15,
            'totalCoveredResidues': 11,
        }
        self.assertEqual(expected, loads(out.getvalue()))

    def testChecksumDirtyOnEmptyDatabase(self):
        """
        The database checksum status must be dirty on creation.
        """
        db = Database([], [])
        self.assertTrue(db.checksumDirty())

    def testChecksumCleanAfterCallingChecksum(self):
        """
        The database checksum status must be non-dirty once the checksum()
        function has been called.
        """
        db = Database([], [])
        db.checksum()
        self.assertFalse(db.checksumDirty())

    def testChecksumCDirtyAfterCallingChecksumThenAddingSubject(self):
        """
        The database checksum status must be dirty if a new subject is added
        after the checksum() function has been called.
        """
        db = Database([], [])
        db.checksum()
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertTrue(db.checksumDirty())

    def testChecksumEmptyDatabase(self):
        """
        The database checksum must be as expected when the database has no
        finders and no subjects.
        """
        db = Database([], [])
        expected = sha256()
        self._update(None, expected)  # Limit per landmark value.
        self._update(None, expected)  # Max distance value.
        self.assertEqual(expected.hexdigest(), db.checksum())

    def testChecksumEmptyDatabaseWithNonDefaultParams(self):
        """
        The database checksum must be as expected when the database is given
        non-default values for limitPerLandmark and maxDistance.
        """
        db = Database([], [], limitPerLandmark=3, maxDistance=10)
        expected = sha256()
        self._update(3, expected)
        self._update(10, expected)
        self.assertEqual(expected.hexdigest(), db.checksum())

    def testChecksumEmptyDatabaseWithFinders(self):
        """
        The database checksum must be as expected when the database has
        finders.
        """
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        expected = sha256()
        self._update(AlphaHelix.NAME, expected)
        self._update(BetaStrand.NAME, expected)
        self._update(Peaks.NAME, expected)
        self._update(Troughs.NAME, expected)
        self._update(None, expected)  # Limit per landmark value.
        self._update(None, expected)  # Max distance value.
        self.assertEqual(expected.hexdigest(), db.checksum())

    def testChecksumEmptyDatabaseWithFinderOrderInvariant(self):
        """
        The database checksum must be identical when the database has finders,
        no matter what order the finders are given.
        """
        db1 = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db2 = Database([BetaStrand, AlphaHelix], [Troughs, Peaks])
        self.assertEqual(db1.checksum(), db2.checksum())

    def testChecksumOneSubjectNoLandmarks(self):
        """
        The database checksum must be as expected when the database has one
        subject but no landmarks are found.
        """
        db = Database([AlphaHelix], [])
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        db.addSubject(subject)

        expected = sha256()
        self._update(AlphaHelix.NAME, expected)
        self._update(None, expected)  # Limit per landmark value.
        self._update(None, expected)  # Max distance value.
        self._update('subject', expected)
        self._update('FRRRFRRRFASAASA', expected)
        self.assertEqual(expected.hexdigest(), db.checksum())

    def testChecksumOneSubjectTwoLandmarks(self):
        """
        The database checksum must be as expected when the database has one
        subject with two landmarks.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        db = Database([AlphaHelix], [])
        subject = AARead('id', sequence)
        db.addSubject(subject)

        expected = sha256()
        self._update(AlphaHelix.NAME, expected)
        self._update(None, expected)  # Limit per landmark value.
        self._update(None, expected)  # Max distance value.
        self._update('id', expected)
        self._update(sequence, expected)

        # The database dictionary now contains
        #
        # {
        #     'A3:A2:31': [{'subjectIndex': 0, 'offset': 31}],
        #     'A2:A3:-31': [{'subjectIndex': 0, 'offset': 0}]
        # }
        #
        # We have to add the keys in order and the keys and values of the
        # values in order.

        self._update('A2:A3:-31', expected)
        self._update('offset', expected)
        self._update(0, expected)
        self._update('subjectIndex', expected)
        self._update(0, expected)

        self._update('A3:A2:31', expected)
        self._update('offset', expected)
        self._update(31, expected)
        self._update('subjectIndex', expected)
        self._update(0, expected)

        self.assertEqual(expected.hexdigest(), db.checksum())
