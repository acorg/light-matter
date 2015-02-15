from unittest import TestCase
from json import loads
from cStringIO import StringIO
from binascii import crc32

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
    def _checksum(strings):
        """
        Compute a database checksum.

        Note that this function knows the details of the database checksum
        algorithm, which is not ideal.

        @param strings: A C{list} of strings to compute the checksum for.
        @return: An C{int} checksum.
        """
        return crc32('\0'.join(map(str, strings)) + '\0', 0x0) & 0xFFFFFFFF

    def testSignificanceFractionDefault(self):
        """
        The significanceFraction default value must be as expected.
        """
        self.assertEqual(0.25, Database.SIGNIFICANCE_FRACTION_DEFAULT)

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

    def testKeyWithFeatureOnLeft(self):
        """
        The database key function must return the expected (negative offset)
        key when the second feature is to the left of the first.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        self.assertEqual('A:B:-10', db.key(landmark, trigPoint))

    def testKeyWithFeatureOnRight(self):
        """
        The database key function must return the expected (positive offset)
        key when the second feature is to the right of the first.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        self.assertEqual('A:B:10', db.key(landmark, trigPoint))

    def testKeyWithSymbolDetail(self):
        """
        The database key function must return the expected value when the
        landmark it is passed has a repeat count.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0, 5)
        trigPoint = TrigPoint('name', 'B', 30)
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
        If one read is added and it has two landmarks, one key is added
        to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:23': {'0': [0]},
            },
            db.d)

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one read is added, the database statistics must be correct.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(36, db.totalResidues)
        self.assertEqual(22, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarksSameOffsets(self):
        """
        If two identical reads are added, both with two landmarks at the same
        offsets, only one key is added to the database and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the database since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:23': {'0': [0], '1': [0]},
            },
            db.d)

    def testTwoReadsTwoLandmarksStatistics(self):
        """
        If two identical reads are added, the database statistics must be
        correct.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(2, db.subjectCount)
        self.assertEqual(72, db.totalResidues)
        self.assertEqual(44, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarksLimitZeroPairsPerLandmark(self):
        """
        If two identical reads are added, both with two landmarks, no keys
        will be added to the dictionary if limitPerLandmark is zero.
        """
        db = Database([AlphaHelix], [], limitPerLandmark=0)
        db.addSubject(AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual({}, db.d)

    def testTwoReadsTwoLandmarksDifferentOffsets(self):
        """
        If two reads are added, both with two landmarks separated by the same
        distance, only one key is added to the database and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the database since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'AFRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2',  'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(
            {
                'A2:A3:23': {'0': [1], '1': [0]},
            },
            db.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one read is added and it has one landmark and one peak, one pair is
        added to the database.
        """
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual(
            {
                'A2:P:10': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkOnePeakBucketFactor(self):
        """
        If a bucketFactor is used, the right distance needs to be calculated.
        The offsets are 10 aa apart, the bucketFactor should make that 2.
        """
        db = Database([AlphaHelix], [Peaks], bucketFactor=5)
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual(
            {
                'A2:P:2': {'0': [0]},
            },
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
        self.assertEqual(
            {
                'A2:P:13': {'0': [0]},
                'A2:P:10': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksLimitOnePairPerLandmark(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        limit of one pair per landmarks is imposed, only one key is added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], limitPerLandmark=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual(
            {
                'A2:P:10': {'0': [0]},
            },
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
        self.assertEqual(
            {
                'A2:P:10': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksLargeMaxDistance(self):
        """
        If one read is added and it has one landmark and two peaks, and a
        maximum distance is imposed that is greater than the distance to the
        peaks, two keys are added to the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=15)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual(
            {
                'A2:P:13': {'0': [0]},
                'A2:P:10': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksPermissiveMinDistance(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        permissive minimum distance is imposed, all keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual(
            {
                'A2:P:13': {'0': [0]},
                'A2:P:10': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMinDistance(self):
        """
        If one read is added and it has one landmark and two peaks, but an
        intermediate minimum distance is imposed, only the key for the pair
        that exceeds the minimum distance is added to the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=11)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual(
            {
                'A2:P:13': {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksSevereMinDistance(self):
        """
        If one read is added and it has one landmark and two peaks, but a
        severe minimum distance is imposed, no keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=100)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({}, db.d)

    def testMultipleSubjectOffsets(self):
        """
        If one read is added and it has one landmark and one peak separated by
        10 bases and then, later in the subject, the same pair with the
        same separation, one key must be added to the database and it
        should have two subject offsets.  Note that minDistance and
        maxDistance are used to discard the matches some longer and shorter
        distance pairs that only have one subject offset (i.e., that only
        appear in the subject once).
        """
        seq = 'FRRRFRRRFASA'
        db = Database([AlphaHelix], [Peaks], minDistance=5, maxDistance=10)
        db.addSubject(AARead('id', seq + seq))
        self.assertEqual(
            {
                'A2:P:10': {'0': [0, 12]},
            },
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
        self.assertEqual(None, result.minDistance)
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
        self.assertEqual(db.minDistance, result.minDistance)
        self.assertEqual(db.trigPointFinderClasses,
                         result.trigPointFinderClasses)
        self.assertEqual(db.totalResidues, result.totalResidues)
        self.assertEqual(db.checksum, result.checksum)

    def testChecksumAfterSaveLoad(self):
        """
        A database that has a sequence added to it, which is then saved and
        then re-loaded, and then has a second sequence is added to it must have
        the same checksum as a database that simply has the two sequences added
        to it without interruption.
        """
        seq1 = 'FRRRFRRRFASAASA'
        seq2 = 'MMMMMMMMMFRRRFR'
        db1 = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db1.addSubject(AARead('id1', seq1))
        fp = StringIO()
        db1.save(fp)
        fp.seek(0)
        db1 = Database.load(fp)
        db1.addSubject(AARead('id2', seq2))

        db2 = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db2.addSubject(AARead('id1', seq1))
        db2.addSubject(AARead('id2', seq2))

        self.assertEqual(db1.checksum, db2.checksum)

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
        default value of significanceFraction.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual(
            {
                0: [
                    {'trigPointName': 'Peaks',
                     'landmarkLength': 9,
                     'readOffset': 0,
                     'subjectLength': 16,
                     'subjectOffsets': [1],
                     'landmarkName': 'AlphaHelix'},
                    {'trigPointName': 'Peaks',
                     'landmarkLength': 9,
                     'readOffset': 0,
                     'subjectLength': 16,
                     'subjectOffsets': [1],
                     'landmarkName': 'AlphaHelix'}
                ]
            }, result.matches)
        self.assertEqual(0, len(list(result.significant())))

    def testFindOneMatchingSignificant(self):
        """
        One matching and significant subject must be found if the
        significanceFraction is sufficiently low.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASA')
        query = AARead('query', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(subject)
        result = db.find(query, significanceFraction=0.0)
        self.assertEqual(
            {
                0: [
                    {
                        'trigPointName': 'Peaks',
                        'landmarkLength': 9,
                        'readOffset': 0,
                        'subjectOffsets': [1],
                        'landmarkName': 'AlphaHelix',
                        'subjectLength': 16
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
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'readOffset': 0,
                        'subjectOffsets': [0],
                        'trigPointName': 'Peaks',
                        'subjectLength': 15,
                    },
                    {
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'readOffset': 0,
                        'subjectOffsets': [0],
                        'trigPointName': 'Peaks',
                        'subjectLength': 15,
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
                      maxDistance=9, minDistance=5)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))

        checksum = self._checksum([
            AlphaHelix.NAME, Peaks.NAME,
            3,  # Limit per landmark.
            9,  # Max distance.
            5,  # Min distance.
            1,  # Bucket factor.
            'id', 'FRRRFRRRFASAASA',
        ])

        out = StringIO()
        db.saveParamsAsJSON(out)
        expected = {
            'checksum': checksum,
            'landmarkFinderClasses': ['AlphaHelix'],
            'trigPointFinderClasses': ['Peaks'],
            'limitPerLandmark': 3,
            'maxDistance': 9,
            'minDistance': 5,
            'subjectCount': 1,
            'totalResidues': 15,
            'totalCoveredResidues': 11,
            'bucketFactor': 1,
        }
        self.assertEqual(expected, loads(out.getvalue()))

    def testChecksumEmptyDatabase(self):
        """
        The database checksum must be as expected when the database has no
        finders and no subjects.
        """
        db = Database([], [])
        checksum = self._checksum([
            None,  # Limit per landmark.
            None,  # Max distance.
            None,  # Min distance.
            1,  # Bucket factor.
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumEmptyDatabaseWithNonDefaultParams(self):
        """
        The database checksum must be as expected when the database is given
        non-default values for limitPerLandmark and maxDistance.
        """
        db = Database([], [], limitPerLandmark=3, maxDistance=9,
                      minDistance=1)
        checksum = self._checksum([
            3,  # Limit per landmark.
            9,  # Max distance.
            1,  # Min distance.
            1,  # Bucket factor.
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumEmptyDatabaseWithFinders(self):
        """
        The database checksum must be as expected when the database has
        finders.
        """
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        checksum = self._checksum([
            AlphaHelix.NAME, BetaStrand.NAME, Peaks.NAME, Troughs.NAME,
            None,  # Limit per landmark.
            None,  # Max distance.
            None,  # Min distance.
            1,  # Bucket factor.
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumEmptyDatabaseWithFinderOrderInvariant(self):
        """
        The database checksum must be identical when the database has finders,
        no matter what order the finders are given.
        """
        db1 = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db2 = Database([BetaStrand, AlphaHelix], [Troughs, Peaks])
        self.assertEqual(db1.checksum, db2.checksum)

    def testChecksumOneSubjectNoLandmarks(self):
        """
        The database checksum must be as expected when the database has one
        subject but no landmarks are found.
        """
        db = Database([AlphaHelix], [])
        subject = AARead('id', 'FRRRFRRRFASAASA')
        db.addSubject(subject)

        checksum = self._checksum([
            AlphaHelix.NAME,
            None,  # Limit per landmark.
            None,  # Max distance.
            None,  # Min distance.
            1,  # Bucket factor.
            'id', 'FRRRFRRRFASAASA',
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumOneSubjectTwoLandmarks(self):
        """
        The database checksum must be as expected when the database has one
        subject with two landmarks.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        db = Database([AlphaHelix], [])
        subject = AARead('id', sequence)
        db.addSubject(subject)
        checksum = self._checksum([
            AlphaHelix.NAME,
            None,  # Limit per landmark.
            None,  # Max distance.
            None,  # Min distance.
            1,  # Bucket factor.
            'id', sequence,
        ])
        self.assertEqual(checksum, db.checksum)

    def testSaveLoadWithNonDefaultParameters(self):
        """
        When asked to save and then load a database with non-default
        parameters, a database with the correct parameters must result.
        """
        db = Database([], [], limitPerLandmark=16, maxDistance=17,
                      minDistance=18, bucketFactor=19)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.load(fp)
        self.assertEqual(16, result.limitPerLandmark)
        self.assertEqual(17, result.maxDistance)
        self.assertEqual(18, result.minDistance)
        self.assertEqual(19, result.bucketFactor)

    def testBucketFactorValueError(self):
        """
        If the bucketFactor is <= 0, a ValueError must be raised.
        """
        error = 'bucketFactor must be > 0\\.'
        self.assertRaisesRegexp(ValueError, error, Database, [], [],
                                bucketFactor=0)
