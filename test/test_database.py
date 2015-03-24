from unittest import TestCase
from json import loads
from cStringIO import StringIO
from binascii import crc32

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark, TrigPoint
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.database import Database, Subject
from light.reads import ScannedRead


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
        self.assertEqual(0.25, Database.DEFAULT_SIGNIFICANCE_FRACTION)

    def testFindersAreStored(self):
        """
        The list of landmark and trig point finders must be stored correctly.
        """
        db = Database([AlphaHelix], [Peaks])
        self.assertEqual([AlphaHelix], db.landmarkClasses)
        self.assertEqual([Peaks], db.trigPointClasses)

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
        distanceMinus10 = str(scale(-10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distanceMinus10,
                         db.hash(landmark, trigPoint))

    def testKeyWithFeatureOnRight(self):
        """
        The database key function must return the expected (positive offset)
        key when the second feature is to the right of the first.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distance10, db.hash(landmark, trigPoint))

    def testKeyWithFeatureOnLeftAndNonDefaultDistanceBase(self):
        """
        The database key function must return the expected key when the
        database has a non-default distance base and the second feature is to
        the left of the first.
        """
        db = Database([], [], distanceBase=1.5)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        distanceMinus10 = str(scale(-10, 1.5))
        self.assertEqual('A:B:' + distanceMinus10,
                         db.hash(landmark, trigPoint))

    def testKeyWithFeatureOnRightAndNonDefaultDistanceBase(self):
        """
        The database key function must return the expected key when the
        database has a non-default distance base and the second feature is to
        the right of the first.
        """
        db = Database([], [], distanceBase=1.5)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, 1.5))
        self.assertEqual('A:B:' + distance10, db.hash(landmark, trigPoint))

    def testKeyWithSymbolDetail(self):
        """
        The database key function must return the expected value when the
        landmark it is passed has a repeat count.
        """
        db = Database([], [])
        landmark = Landmark('name', 'A', 20, 0, 5)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A5:B:' + distance10, db.hash(landmark, trigPoint))

    def testInitialDatabaseHasNoSubjectInfo(self):
        """
        The database must not have any stored subject information if no
        subjects have been added.
        """
        db = Database([], [])
        self.assertEqual([], list(db.getSubjects()))

    def testAddSubjectReturnsIndex(self):
        """
        If one subject is added, addSubject must return the index (0) of the
        added subject.
        """
        db = Database([AlphaHelix], [])
        self.assertEqual(0, db.addSubject(AARead('id', 'FRRRFRRRF')))

    def testAddSameSubjectReturnsSameIndex(self):
        """
        If an identical subject is added multiple times, the same subject
        index must be returned.
        """
        db = Database([AlphaHelix], [])
        self.assertEqual(db.addSubject(AARead('id', 'FRRRFRRRF')),
                         db.addSubject(AARead('id', 'FRRRFRRRF')))

    def testAddSameSubjectLeavesDatabaseSizeTheSame(self):
        """
        If an identical subject is added multiple times, the database size
        does not increase.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)

    def testGetSubjectIndexError(self):
        """
        If an out-of-range subject index is passed to getSubject, an IndexError
        must be raised.
        """
        db = Database([AlphaHelix], [])
        error = '^list index out of range$'
        self.assertRaisesRegexp(IndexError, error, db.getSubject, 0)
        self.assertRaisesRegexp(IndexError, error, db.getSubject, -1)

    def testGetSubjectKeyError(self):
        """
        If a non-existent subject is passed to getSubject, a KeyError
        must be raised.
        """
        db = Database([AlphaHelix], [])
        error = "^'id'$"
        self.assertRaisesRegexp(KeyError, error, db.getSubject,
                                AARead('id', 'FF'))

    def testGetSubjectByIndex(self):
        """
        If a subject is added, getSubject must be able to return it given an
        integer index.
        """
        db = Database([AlphaHelix], [])
        subject = AARead('id', 'FRRRFRRRF')
        index = db.addSubject(subject)
        self.assertEqual(AARead('id', 'FRRRFRRRF'), db.getSubject(index))

    def testGetSubjectBySubject(self):
        """
        If a subject is added, getSubject must be able to return it given an
        identical subject to look up.
        """
        db = Database([AlphaHelix], [])
        index = db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(index, db.getSubject(AARead('id', 'FRRRFRRRF')))

    def testGetSubjectHashCount(self):
        """
        If a subject is added, getSubject must return a Subject instance that
        has the correct hash count.
        """
        db = Database([AlphaHelix], [])
        subject = AARead('id', 'FRRRFRRRFAFRRRFRRRF')
        index = db.addSubject(subject)
        self.assertEqual(1, db.getSubject(index).hashCount)

    def testOneReadOneLandmark(self):
        """
        If one subject is added but it only has one landmark, nothing is added
        to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkGetSubjects(self):
        """
        If one subject with just one landmark (and hence no hashes) is added,
        an entry is appended to the database subject info.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        subjects = list(db.getSubjects())
        self.assertEqual(1, len(subjects))
        subject = subjects[0]
        self.assertEqual(AARead('id', 'FRRRFRRRF'), subject)
        self.assertEqual(0, subject.hashCount)

    def testOneReadTwoLandmarksGetSubjects(self):
        """
        If one subject with two landmarks (and hence one hash) is added, an
        entry is appended to the database subject info.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFAFRRRFRRRF'))
        subject = list(db.getSubjects())[0]
        self.assertEqual(AARead('id', 'FRRRFRRRFAFRRRFRRRF'), subject)
        self.assertEqual(1, subject.hashCount)

    def testOneReadOneLandmarkStatistics(self):
        """
        If one subject is added the database statistics must be correct.
        """
        db = Database([], [])
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(9, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testOneReadTwoLandmarks(self):
        """
        If one subject is added and it has two landmarks, one key is added
        to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        distance23 = str(scale(23, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [0]},
            },
            db.d)

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one subject is added, the database statistics must be correct.
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
        distance23 = str(scale(23, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [0], '1': [0]},
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
        If two subjects are added, both with two landmarks separated by the
        same distance, only one key is added to the database and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the database since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id1', 'AFRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        distance23 = str(scale(23, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [1], '1': [0]},
            },
            db.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one subject is added and it has one landmark and one peak, one pair
        is added to the database.
        """
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkOnePeakDistanceBase(self):
        """
        If a non-default distanceBase of 2.0 is used, the right distance needs
        to be calculated. In this case, the offsets are 10 AA apart, and the
        distanceBase scaling will change that to a 3 (since int(log base 2 10)
        = 3), though we don't test the 3 value explicitly since that may change
        if we ever change the scale function. That's desirable, but we already
        have tests in test_distance.py that will break in that case.
        """
        distanceBase = 2.0
        db = Database([AlphaHelix], [Peaks], distanceBase=distanceBase)
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        distance10 = str(scale(10, distanceBase))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkOnePeakNoTrigFinders(self):
        """
        If one subject is added and it has one landmark and one peak, but no
        trig finders are in use, nothing is added to the database.
        """
        db = Database([AlphaHelix], [])
        db.addSubject(AARead('id', 'FRRRFRRRFASA'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkTwoPeaks(self):
        """
        If one subject is added and it has one landmark and two peaks, two
        pairs are added to the database.
        """
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance13 = str(scale(13, Database.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksLimitOnePairPerLandmark(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        limit of one pair per landmarks is imposed, only one key is added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], limitPerLandmark=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksSevereMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe maximum distance is imposed, no keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({}, db.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        maximum distance is imposed that makes one of the peaks too far
        away, only one key is added to the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksLargeMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, and a
        maximum distance is imposed that is greater than the distance to the
        peaks, two keys are added to the database.
        """
        db = Database([AlphaHelix], [Peaks], maxDistance=15)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance13 = str(scale(13, Database.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksPermissiveMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        permissive minimum distance is imposed, all keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=1)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance13 = str(scale(13, Database.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but an
        intermediate minimum distance is imposed, only the key for the pair
        that exceeds the minimum distance is added to the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=11)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        distance13 = str(scale(13, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
            },
            db.d)

    def testOneReadOneLandmarkTwoPeaksSevereMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe minimum distance is imposed, no keys are added to
        the database.
        """
        db = Database([AlphaHelix], [Peaks], minDistance=100)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        self.assertEqual({}, db.d)

    def testMultipleSubjectOffsets(self):
        """
        If one subject is added and it has one landmark and one peak separated
        by 10 bases and then, later in the subject, the same pair with the
        same separation, one key must be added to the database and it
        should have two subject offsets.  Note that minDistance and
        maxDistance are used to discard the matches some longer and shorter
        distance pairs that only have one subject offset (i.e., that only
        appear in the subject once).
        """
        seq = 'FRRRFRRRFASA'
        db = Database([AlphaHelix], [Peaks], minDistance=5, maxDistance=10)
        db.addSubject(AARead('id', seq + seq))
        distance10 = str(scale(10, Database.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0, 12]},
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
        self.assertEqual([], result.landmarkClasses)
        self.assertEqual(Database.DEFAULT_LIMIT_PER_LANDMARK,
                         result.limitPerLandmark)
        self.assertEqual(Database.DEFAULT_MAX_DISTANCE, result.maxDistance)
        self.assertEqual(Database.DEFAULT_MIN_DISTANCE, result.minDistance)
        self.assertEqual([], result.trigPointClasses)
        self.assertEqual(0, result.totalResidues)

    def testSaveLoadMissingLandmark(self):
        """
        If a database is saved with a landmark class whose name cannot be
        found when the database is later loaded, a ValueError error must be
        raised.
        """
        db = Database([AlphaHelix], [])
        fp = StringIO()
        db.save(fp)
        newSave = fp.getvalue().replace('AlphaHelix', 'Non-existent')
        error = ('^Could not find landscape finder class Non-existent. '
                 'Has that class been renamed or removed\?$')
        self.assertRaisesRegexp(ValueError, error, db.load, StringIO(newSave))

    def testSaveLoadMissingTrigPoint(self):
        """
        If a database is saved with a trig point class whose name cannot be
        found when the database is later loaded, a ValueError error must be
        raised.
        """
        db = Database([], [Peaks])
        fp = StringIO()
        db.save(fp)
        newSave = fp.getvalue().replace('Peaks', 'Non-existent')
        error = ('^Could not find trig point finder class Non-existent. '
                 'Has that class been renamed or removed\?$')
        self.assertRaisesRegexp(ValueError, error, db.load, StringIO(newSave))

    def testLoadInvalidJSON(self):
        """
        If a database is attempted from a file that does not contain valid
        JSON, a ValueError error must be raised.
        """
        db = Database([], [Peaks])
        error = '^Expected object or value$'
        self.assertRaisesRegexp(ValueError, error, db.load, StringIO('xxx'))

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
        self.assertEqual(db._subjectInfo, result._subjectInfo)
        self.assertEqual(db._idSequenceCache, result._idSequenceCache)
        self.assertEqual(db.landmarkClasses, result.landmarkClasses)
        self.assertEqual(db.limitPerLandmark, result.limitPerLandmark)
        self.assertEqual(db.maxDistance, result.maxDistance)
        self.assertEqual(db.minDistance, result.minDistance)
        self.assertEqual(db.trigPointClasses, result.trigPointClasses)
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
        subject = AARead('subject', 'AFRRRFRRRFASAASAVVVVVVASAVVVASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')
        db = Database([AlphaHelix, BetaStrand], [Peaks])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual(
            {
                0: [
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryLandmarkOffsets': [0],
                        'queryTrigPointOffsets': [10],
                        'subjectLandmarkOffsets': [1],
                        'trigPoint': TrigPoint('Peaks', 'P', 10),
                    },
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryLandmarkOffsets': [0],
                        'queryTrigPointOffsets': [13],
                        'subjectLandmarkOffsets': [1],
                        'trigPoint': TrigPoint('Peaks', 'P', 13),
                    }
                ]
            }, result.matches)
        self.assertEqual(0, len(list(result.significant())))

    def testFindOneMatchingSignificant(self):
        """
        One matching and significant subject must be found if the
        significanceFraction is sufficiently low.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        db = Database([AlphaHelix], [Peaks], maxDistance=11)
        db.addSubject(subject)
        result = db.find(query, significanceFraction=0.0)
        self.assertEqual(
            {
                0: [
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                        'queryLandmarkOffsets': [1],
                        'queryTrigPointOffsets': [11],
                        'subjectLandmarkOffsets': [1],
                        'trigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                ],
            },
            result.matches)

    def testSymmetricFindScoresSameSubjectAndQuery(self):
        """
        The score of matching a sequence A against a sequence B must
        be the same as when matching B against A, and that score must
        be 1.0 when the subject and the query are identical.
        """
        sequence = 'AFRRRFRRRFASAASAFRRRFRRRF'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        db = Database([AlphaHelix, BetaStrand], [Peaks])
        db.addSubject(subject)
        result = db.find(query, significanceFraction=0.0)
        score1 = result.analysis[0]['bestScore']

        db = Database([AlphaHelix, BetaStrand], [Peaks])
        db.addSubject(query)
        result = db.find(subject, significanceFraction=0.0)
        score2 = result.analysis[0]['bestScore']

        self.assertEqual(score1, score2)
        self.assertEqual(1.0, score1)

    def testSymmetricFindScoresDifferingSubjectAndQuery(self):
        """
        The score of matching a sequence A against a sequence B must
        be the same as when matching B against A, including when the number
        of hashes in the two differs and the scores are not 1.0.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASAFRRRFRRRF')
        query = AARead('query', 'FRRRFRRRFASAVVVVVV')
        db = Database([AlphaHelix, BetaStrand], [Peaks])
        db.addSubject(subject)
        hashCount1 = db.getSubject(0).hashCount
        result = db.find(query, significanceFraction=0.0)
        score1 = result.analysis[0]['bestScore']

        db = Database([AlphaHelix, BetaStrand], [Peaks])
        db.addSubject(query)
        hashCount2 = db.getSubject(0).hashCount
        result = db.find(subject, significanceFraction=0.0)
        score2 = result.analysis[0]['bestScore']

        self.assertNotEqual(hashCount1, hashCount2)
        self.assertEqual(score1, score2)
        self.assertNotEqual(1.0, score1)

    def testFindNoneMatchingTooSmallDistance(self):
        """
        One matching key must be found.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        db = Database([AlphaHelix], [Peaks], maxDistance=1)
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindNoneMatchingNoTrigPoint(self):
        """
        One matching key must be found.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        db = Database([AlphaHelix], [])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindTwoMatchingInSameSubject(self):
        """
        Two matching hashes in the subject must be found correctly.
        """
        sequence = 'FRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual(
            {
                0: [
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryLandmarkOffsets': [0],
                        'queryTrigPointOffsets': [10],
                        'subjectLandmarkOffsets': [0],
                        'trigPoint': TrigPoint('Peaks', 'P', 10),
                    },
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryLandmarkOffsets': [0],
                        'queryTrigPointOffsets': [13],
                        'subjectLandmarkOffsets': [0],
                        'trigPoint': TrigPoint('Peaks', 'P', 13),
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
                      maxDistance=19, minDistance=5)
        sequence = 'AFRRRFRRRFASAASA'
        db.addSubject(AARead('id', sequence))

        checksum = self._checksum([
            AlphaHelix.NAME,
            AlphaHelix.SYMBOL,
            Peaks.NAME,
            Peaks.SYMBOL,
            3,  # Limit per landmark.
            19,  # Max distance.
            5,  # Min distance.
            Database.DEFAULT_DISTANCE_BASE,
            'id',
            sequence,
            2,  # Hash count.
        ])

        out = StringIO()
        db.saveParamsAsJSON(out)
        expected = {
            'checksum': checksum,
            'landmarkClasses': ['AlphaHelix'],
            'trigPointClasses': ['Peaks'],
            'limitPerLandmark': 3,
            'maxDistance': 19,
            'minDistance': 5,
            'subjectCount': 1,
            'totalResidues': 16,
            'totalCoveredResidues': 11,
            'distanceBase': Database.DEFAULT_DISTANCE_BASE,
        }
        self.assertEqual(expected, loads(out.getvalue()))

    def testChecksumEmptyDatabase(self):
        """
        The database checksum must be as expected when the database has no
        finders and no subjects.
        """
        db = Database([], [])
        checksum = self._checksum([
            Database.DEFAULT_LIMIT_PER_LANDMARK,
            Database.DEFAULT_MAX_DISTANCE,
            Database.DEFAULT_MIN_DISTANCE,
            Database.DEFAULT_DISTANCE_BASE,
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumEmptyDatabaseWithNonDefaultParams(self):
        """
        The database checksum must be as expected when the database is given
        non-default values for limitPerLandmark and maxDistance.
        """
        db = Database([], [], limitPerLandmark=3, maxDistance=9,
                      minDistance=1, distanceBase=1.5)
        checksum = self._checksum([
            3,  # Limit per landmark.
            9,  # Max distance.
            1,  # Min distance.
            1.5,  # Distance base.
        ])
        self.assertEqual(checksum, db.checksum)

    def testChecksumEmptyDatabaseWithFinders(self):
        """
        The database checksum must be as expected when the database has
        finders.
        """
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs])
        checksum = self._checksum([
            AlphaHelix.NAME,
            BetaStrand.NAME,
            AlphaHelix.SYMBOL,
            BetaStrand.SYMBOL,
            Peaks.NAME,
            Troughs.NAME,
            Peaks.SYMBOL,
            Troughs.SYMBOL,
            Database.DEFAULT_LIMIT_PER_LANDMARK,
            Database.DEFAULT_MAX_DISTANCE,
            Database.DEFAULT_MIN_DISTANCE,
            Database.DEFAULT_DISTANCE_BASE,
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
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        db.addSubject(subject)

        checksum = self._checksum([
            AlphaHelix.NAME,
            AlphaHelix.SYMBOL,
            Database.DEFAULT_LIMIT_PER_LANDMARK,
            Database.DEFAULT_MAX_DISTANCE,
            Database.DEFAULT_MIN_DISTANCE,
            Database.DEFAULT_DISTANCE_BASE,
            'id',
            sequence,
            0,  # Hash count.
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
            AlphaHelix.SYMBOL,
            Database.DEFAULT_LIMIT_PER_LANDMARK,
            Database.DEFAULT_MAX_DISTANCE,
            Database.DEFAULT_MIN_DISTANCE,
            Database.DEFAULT_DISTANCE_BASE,
            'id',
            sequence,
            1,  # Hash count.
        ])
        self.assertEqual(checksum, db.checksum)

    def testSaveLoadWithNonDefaultParameters(self):
        """
        When asked to save and then load a database with non-default
        parameters, a database with the correct parameters must result.
        """
        db = Database([], [], limitPerLandmark=16, maxDistance=17,
                      minDistance=18, distanceBase=19.0)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.load(fp)
        self.assertEqual(16, result.limitPerLandmark)
        self.assertEqual(17, result.maxDistance)
        self.assertEqual(18, result.minDistance)
        self.assertEqual(19.0, result.distanceBase)

    def testDistanceBaseZeroValueError(self):
        """
        If the distanceBase is zero, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        self.assertRaisesRegexp(ValueError, error, Database, [], [],
                                distanceBase=0.0)

    def testDistanceBaseLessThanZeroValueError(self):
        """
        If the distanceBase is < 0, a ValueError must be raised.
        """
        error = 'distanceBase must be > 0\\.'
        self.assertRaisesRegexp(ValueError, error, Database, [], [],
                                distanceBase=-1.0)

    def testScan(self):
        """
        The scan method must return a scanned subject.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks])
        db.addSubject(subject)
        scannedSubject = db.scan(subject)
        self.assertIsInstance(scannedSubject, ScannedRead)

    def testGetScannedPairs(self):
        """
        The getSequencePairs method must return pairs of
        (landmark, trigPoints).
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix], [Peaks], distanceBase=1.0)
        db.addSubject(subject)
        scannedSubject = db.scan(subject)
        pairs = list(db.getScannedPairs(scannedSubject))
        # First pair.
        landmark, trigPoint = pairs[0]
        self.assertEqual(Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL,
                                  0, 9, 2), landmark)
        self.assertEqual(TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10), trigPoint)
        # Second pair.
        landmark, trigPoint = pairs[1]
        self.assertEqual(Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL,
                                  0, 9, 2), landmark)
        self.assertEqual(TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13), trigPoint)
        self.assertEqual(2, len(pairs))

    def testCollectReadHashesWithNoHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets. The result should be
        empty if there are no landmarks in the read.
        """
        db = Database([AlphaHelix], [])
        query = AARead('query', 'AAA')
        scannedQuery = db.scan(query)
        hashCount = db.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashesWithOneLandmark(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets. The result should be
        empty if there is only one landmark in the read.
        """
        db = Database([AlphaHelix], [])
        query = AARead('query', 'FRRRFRRRF')
        scannedQuery = db.scan(query)
        hashCount = db.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets.
        """
        db = Database([AlphaHelix], [Peaks], distanceBase=1.0)
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFASAASA')
        scannedQuery = db.scan(query)
        hashCount = db.getHashes(scannedQuery)
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        helixAt15 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 15, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)
        peakAt13 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13)
        peakAt25 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 25)
        peakAt28 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 28)
        self.assertEqual(
            {
                'A2:A2:15': {
                    'landmark': helixAt0,
                    'landmarkOffsets': [0],
                    'trigPointOffsets': [15],
                    'trigPoint': helixAt15,
                },
                'A2:P:-2': {
                    'landmark': helixAt15,
                    'landmarkOffsets': [15],
                    'trigPointOffsets': [13],
                    'trigPoint': peakAt13,
                },
                'A2:P:-5': {
                    'landmark': helixAt15,
                    'landmarkOffsets': [15],
                    'trigPointOffsets': [10],
                    'trigPoint': peakAt10,
                },
                'A2:P:10': {
                    'landmark': helixAt0,
                    'landmarkOffsets': [0, 15],
                    'trigPointOffsets': [10, 25],
                    'trigPoint': peakAt10,
                },
                'A2:P:13': {
                    'landmark': helixAt0,
                    'landmarkOffsets': [0, 15],
                    'trigPointOffsets': [13, 28],
                    'trigPoint': peakAt13,
                },
                'A2:P:25': {
                    'landmark': helixAt0,
                    'landmarkOffsets': [0],
                    'trigPointOffsets': [25],
                    'trigPoint': peakAt25,
                },
                'A2:P:28': {
                    'landmark': helixAt0,
                    'landmarkOffsets': [0],
                    'trigPointOffsets': [28],
                    'trigPoint': peakAt28,
                },
            }, hashCount)

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        fp = StringIO()
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs],
                      limitPerLandmark=16, maxDistance=10, minDistance=0,
                      distanceBase=1)
        db.addSubject(subject)
        db.print_(fp)
        expected = (
            'Landmark finders:\n'
            '  AlphaHelix\n'
            '  BetaStrand\n'
            'Trig point finders:\n'
            '  Peaks\n'
            '  Troughs\n'
            'Subject count: 1\n'
            'Hash count: 3\n'
            'Total residues: 15\n'
            'Coverage: 73.33%\n'
            'Checksum: 682086972\n')
        self.assertEqual(expected, fp.getvalue())

    def testPrintWithHashes(self):
        """
        The print_ function should produce the expected output when asked to
        print hash information.
        """
        fp = StringIO()
        subject = AARead('id', 'FRRRFRRRFASAASA')
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs],
                      limitPerLandmark=16, maxDistance=10, minDistance=0,
                      distanceBase=1)
        db.addSubject(subject)
        db.print_(fp, printHashes=True)
        expected = (
            'Landmark finders:\n'
            '  AlphaHelix\n'
            '  BetaStrand\n'
            'Trig point finders:\n'
            '  Peaks\n'
            '  Troughs\n'
            'Subject count: 1\n'
            'Hash count: 3\n'
            'Total residues: 15\n'
            'Coverage: 73.33%\n'
            'Checksum: 3279991028\n'
            'Subjects (with offsets) by hash:\n'
            '   A2:P:10\n'
            '    id [0]\n'
            '   A2:T:8\n'
            '    id [0]\n'
            '   A2:T:4\n'
            '    id [0]\n'
            'Landmark symbol counts:\n'
            '  AlphaHelix (A2): 3\n'
            'Trig point symbol counts:\n'
            '  Peaks (P): 1\n'
            '  Troughs (T): 2\n')
        self.assertEqual(expected, fp.getvalue())

    def testPrintNoHashes(self):
        """
        The print_ function should report the expected result if no hashes are
        found in the subject.
        """
        fp = StringIO()
        subject = AARead('subject', '')
        db = Database([AlphaHelix, BetaStrand], [Peaks, Troughs],
                      limitPerLandmark=16, maxDistance=10, minDistance=0,
                      distanceBase=1)
        db.addSubject(subject)
        db.print_(fp)
        expected = (
            'Landmark finders:\n'
            '  AlphaHelix\n'
            '  BetaStrand\n'
            'Trig point finders:\n'
            '  Peaks\n'
            '  Troughs\n'
            'Subject count: 1\n'
            'Hash count: 0\n'
            'Total residues: 0\n'
            'Coverage: 0.00%\n'
            'Checksum: 4224788348\n')
        self.assertEqual(expected, fp.getvalue())


class TestSubject(TestCase):
    """
    Tests for the light.database.Subject class.
    """
    def testIsAARead(self):
        """
        A Subject is a subclass of AARead.
        """
        self.assertTrue(isinstance(Subject('id', 'AA', 0), AARead))

    def testAAProperties(self):
        """
        A Subject must call AARead.__init__ with the correct arguments.
        """
        self.assertEqual(AARead('id', 'A', '!'), Subject('id', 'A', 6, '!'))

    def testHashCountIsStored(self):
        """
        A Subject must store the hashcount it is passed.
        """
        self.assertEqual(6, Subject('id', 'AA', 6).hashCount)
