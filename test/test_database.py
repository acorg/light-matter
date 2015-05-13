from unittest import TestCase
from io import StringIO

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark, TrigPoint
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.checksum import Checksum
from light.database import Parameters, Database, Backend, SimpleConnector
from light.reads import ScannedRead


class TestDatabase(TestCase):
    """
    Tests for the light.database.Database class.
    """
    def testSignificanceFractionDefault(self):
        """
        The significanceFraction default value must be as expected.
        """
        self.assertEqual(0.25, Parameters.DEFAULT_SIGNIFICANCE_FRACTION)

    def testFindersAreStored(self):
        """
        The list of landmark and trig point finders must be stored correctly.
        """
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
        self.assertEqual([AlphaHelix], db.params.landmarkClasses)
        self.assertEqual([Peaks], db.params.trigPointClasses)

    def testInitialStatistics(self):
        """
        The database statistics must be initially correct.
        """
        params = Parameters([], [])
        db = Database(params)
        self.assertEqual(0, db.subjectCount)
        self.assertEqual(0, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testKeyWithFeatureOnLeft(self):
        """
        The database key function must return the expected (negative offset)
        key when the second feature is to the left of the first.
        """
        params = Parameters([], [])
        db = Database(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        distanceMinus10 = str(scale(-10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distanceMinus10,
                         db.hash(landmark, trigPoint))

    def testKeyWithFeatureOnRight(self):
        """
        The database key function must return the expected (positive offset)
        key when the second feature is to the right of the first.
        """
        params = Parameters([], [])
        db = Database(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distance10, db.hash(landmark, trigPoint))

    def testKeyWithFeatureOnLeftAndNonDefaultDistanceBase(self):
        """
        The database key function must return the expected key when the
        database has a non-default distance base and the second feature is to
        the left of the first.
        """
        params = Parameters([], [], distanceBase=1.5)
        db = Database(params)
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
        params = Parameters([], [], distanceBase=1.5)
        db = Database(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, 1.5))
        self.assertEqual('A:B:' + distance10, db.hash(landmark, trigPoint))

    def testKeyWithSymbolDetail(self):
        """
        The database key function must return the expected value when the
        landmark it is passed has a repeat count.
        """
        params = Parameters([], [])
        db = Database(params)
        landmark = Landmark('name', 'A', 20, 0, 5)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A5:B:' + distance10, db.hash(landmark, trigPoint))

    def testInitialDatabaseHasNoSubjectInfo(self):
        """
        The database must not have any stored subject information if no
        subjects have been added.
        """
        params = Parameters([], [])
        db = Database(params)
        self.assertEqual([], list(db.getSubjects()))

    def testAddSubjectReturnsIndex(self):
        """
        If one subject is added, addSubject must return the index (0) of the
        added subject.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        self.assertEqual(0, db.addSubject(AARead('id', 'FRRRFRRRF')))

    def testAddSameSubjectReturnsSameIndex(self):
        """
        If an identical subject is added multiple times, the same subject
        index must be returned.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        self.assertEqual(db.addSubject(AARead('id', 'FRRRFRRRF')),
                         db.addSubject(AARead('id', 'FRRRFRRRF')))

    def testAddSameSubjectLeavesDatabaseSizeTheSame(self):
        """
        If an identical subject is added multiple times, the database size
        does not increase.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)

    def testGetSubjectIndexError(self):
        """
        If an out-of-range subject index is passed to getSubject, an IndexError
        must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        error = '^list index out of range$'
        self.assertRaisesRegex(IndexError, error, db.getSubject, 0)
        self.assertRaisesRegex(IndexError, error, db.getSubject, -1)

    def testGetSubjectKeyError(self):
        """
        If a non-existent subject is passed to getSubject, a KeyError
        must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        error = "^'id'$"
        self.assertRaisesRegex(KeyError, error, db.getSubject,
                               AARead('id', 'FF'))

    def testGetSubjectByIndex(self):
        """
        If a subject is added, getSubject must be able to return it given an
        integer index.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        subject = AARead('id', 'FRRRFRRRF')
        index = db.addSubject(subject)
        self.assertEqual(AARead('id', 'FRRRFRRRF'), db.getSubject(index))

    def testGetSubjectBySubject(self):
        """
        If a subject is added, getSubject must be able to return it given an
        identical subject to look up.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        index = db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(index, db.getSubject(AARead('id', 'FRRRFRRRF')))

    def testGetSubjectHashCount(self):
        """
        If a subject is added, getSubject must return a Subject instance that
        has the correct hash count.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        subject = AARead('id', 'FRRRFRRRFAFRRRFRRRF')
        index = db.addSubject(subject)
        self.assertEqual(1, db.getSubject(index).hashCount)

    def testOneReadOneLandmarkGetSubjects(self):
        """
        If one subject with just one landmark (and hence no hashes) is added,
        an entry is appended to the database subject info.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
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
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRFAFRRRFRRRF'))
        subject = list(db.getSubjects())[0]
        self.assertEqual(AARead('id', 'FRRRFRRRFAFRRRFRRRF'), subject)
        self.assertEqual(1, subject.hashCount)

    def testOneReadOneLandmarkStatistics(self):
        """
        If one subject is added the database statistics must be correct.
        """
        params = Parameters([], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(9, db.totalResidues)
        self.assertEqual(0, db.totalCoveredResidues)

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one subject is added, the database statistics must be correct.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(1, db.subjectCount)
        self.assertEqual(36, db.totalResidues)
        self.assertEqual(22, db.totalCoveredResidues)

    def testTwoReadsTwoLandmarksStatistics(self):
        """
        If two identical reads are added, the database statistics must be
        correct.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(2, db.subjectCount)
        self.assertEqual(72, db.totalResidues)
        self.assertEqual(44, db.totalCoveredResidues)

    def testSaveRestoreEmpty(self):
        """
        When asked to save and then restore an empty database, the correct
        database must result.
        """
        params = Parameters([], [])
        db = Database(params)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.restore(fp)
        self.assertEqual(0, result.subjectCount)
        self.assertEqual(0, result.totalCoveredResidues)
        self.assertEqual(0, result.totalResidues)
        self.assertEqual([], result.params.landmarkClasses)
        self.assertEqual(Parameters.DEFAULT_LIMIT_PER_LANDMARK,
                         result.params.limitPerLandmark)
        self.assertEqual(Parameters.DEFAULT_MAX_DISTANCE,
                         result.params.maxDistance)
        self.assertEqual(Parameters.DEFAULT_MIN_DISTANCE,
                         result.params.minDistance)
        self.assertEqual([], result.params.trigPointClasses)

    def testRestoreInvalidJSON(self):
        """
        If a database restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        params = Parameters([], [Peaks])
        db = Database(params)
        error = '^Expected object or value$'
        self.assertRaisesRegex(ValueError, error, db.restore, StringIO('xxx'))

    def testSaveLoadNonEmpty(self):
        """
        When asked to save and then load a non-empty database, the correct
        database must result.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.restore(fp)
        self.assertEqual(db.subjectCount, result.subjectCount)
        self.assertEqual(db.totalCoveredResidues, result.totalCoveredResidues)
        self.assertEqual(db._subjectInfo, result._subjectInfo)
        self.assertEqual(db._idSequenceCache, result._idSequenceCache)
        self.assertEqual(db.params.landmarkClasses,
                         result.params.landmarkClasses)
        self.assertEqual(db.params.limitPerLandmark,
                         result.params.limitPerLandmark)
        self.assertEqual(db.params.maxDistance, result.params.maxDistance)
        self.assertEqual(db.params.minDistance, result.params.minDistance)
        self.assertEqual(db.params.trigPointClasses,
                         result.params.trigPointClasses)
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
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db1 = Database(params1)
        db1.addSubject(AARead('id1', seq1))
        fp = StringIO()
        db1.save(fp)
        fp.seek(0)
        db1 = Database.restore(fp)
        db1.addSubject(AARead('id2', seq2))

        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db2 = Database(params2)
        db2.addSubject(AARead('id1', seq1))
        db2.addSubject(AARead('id2', seq2))

        self.assertEqual(db1.checksum, db2.checksum)

    def testFindNoMatching(self):
        """
        A non-matching key must not be found.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRR')
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
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
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params)
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
        self.assertEqual(0, len(list(result.significantSubjects())))

    def testFindOneMatchingSignificant(self):
        """
        One matching and significant subject must be found if the
        significanceFraction is sufficiently low.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        db = Database(params)
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

    def testFindNoneMatchingTooSmallDistance(self):
        """
        No matches should be found if the max distance is too small.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=1)
        db = Database(params)
        db.addSubject(subject)
        result = db.find(query)
        self.assertEqual({}, result.matches)

    def testFindNoneMatchingNoTrigPoint(self):
        """
        No matches should be found if there is only one landmark and there are
        no trig point finders.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [])
        db = Database(params)
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
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
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

    def testSymmetricFindScoresSameSubjectAndQuery(self):
        """
        The score of matching a sequence A against a sequence B must
        be the same as when matching B against A, and that score must
        be 1.0 when the subject and the query are identical.
        """
        sequence = 'AFRRRFRRRFASAASAFRRRFRRRF'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params)
        db.addSubject(subject)
        result = db.find(query, significanceFraction=0.0)
        score1 = result.analysis[0]['bestScore']

        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params)
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
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params1)
        db.addSubject(subject)
        hashCount1 = db.getSubject(0).hashCount
        result = db.find(query, significanceFraction=0.0)
        score1 = result.analysis[0]['bestScore']

        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params2)
        db.addSubject(query)
        hashCount2 = db.getSubject(0).hashCount
        result = db.find(subject, significanceFraction=0.0)
        score2 = result.analysis[0]['bestScore']

        self.assertNotEqual(hashCount1, hashCount2)
        self.assertEqual(score1, score2)
        self.assertNotEqual(1.0, score1)

    def testChecksumEmptyDatabase(self):
        """
        The database checksum must be the same as the checksum for its
        parameters when no subjects have been added to the database.
        """
        params = Parameters([], [])
        db = Database(params)
        self.assertEqual(params.checksum, db.checksum)

    def testChecksumAfterSubjectAdded(self):
        """
        The database checksum must be as expected when a subject has been
        added to the database.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        db.addSubject(subject)

        checksum = Checksum(params.checksum).update([
            'id',
            sequence,
            '0',  # Hash count.
        ]).checksum
        self.assertEqual(checksum, db.checksum)

    def testSaveLoadWithNonDefaultParameters(self):
        """
        When asked to save and then load a database with non-default
        parameters, a database with the correct parameters must result.
        """
        params = Parameters([], [], limitPerLandmark=16, maxDistance=17,
                            minDistance=18, distanceBase=19.0)
        db = Database(params)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.restore(fp)
        self.assertIs(None, params.compare(result.params))

    def testScan(self):
        """
        The scan method must return a scanned subject.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
        db.addSubject(subject)
        scannedSubject = db.scan(subject)
        self.assertIsInstance(scannedSubject, ScannedRead)

    def testGetScannedPairs(self):
        """
        The getSequencePairs method must return pairs of
        (landmark, trigPoints).
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        db = Database(params)
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
        params = Parameters([AlphaHelix], [])
        db = Database(params)
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
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        query = AARead('query', 'FRRRFRRRF')
        scannedQuery = db.scan(query)
        hashCount = db.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets.
        """
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        db = Database(params)
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
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
        db.print_(fp)
        expected = (
            'Parameters:\n'
            '  Landmark finders:\n'
            '    AlphaHelix\n'
            '    BetaStrand\n'
            '  Trig point finders:\n'
            '    Peaks\n'
            '    Troughs\n'
            '  Limit per landmark: 16\n'
            '  Max distance: 10\n'
            '  Min distance: 0\n'
            '  Distance base: 1\n'
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
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
        db.print_(fp, printHashes=True)
        expected = (
            'Parameters:\n'
            '  Landmark finders:\n'
            '    AlphaHelix\n'
            '    BetaStrand\n'
            '  Trig point finders:\n'
            '    Peaks\n'
            '    Troughs\n'
            '  Limit per landmark: 16\n'
            '  Max distance: 10\n'
            '  Min distance: 0\n'
            '  Distance base: 1\n'
            'Subject count: 1\n'
            'Hash count: 3\n'
            'Total residues: 15\n'
            'Coverage: 73.33%\n'
            'Checksum: 1144016651\n'
            'Backend \'localhost\':\n'
            'Hash count: 3\n'
            'Checksum: 1144016651\n'
            'Subjects (with offsets) by hash:\n'
            '   A2:P:10\n'
            '    0 [0]\n'
            '   A2:T:4\n'
            '    0 [0]\n'
            '   A2:T:8\n'
            '    0 [0]\n'
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
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
        db.print_(fp)
        expected = (
            'Parameters:\n'
            '  Landmark finders:\n'
            '    AlphaHelix\n'
            '    BetaStrand\n'
            '  Trig point finders:\n'
            '    Peaks\n'
            '    Troughs\n'
            '  Limit per landmark: 16\n'
            '  Max distance: 10\n'
            '  Min distance: 0\n'
            '  Distance base: 1\n'
            'Subject count: 1\n'
            'Hash count: 0\n'
            'Total residues: 0\n'
            'Coverage: 0.00%\n'
            'Checksum: 4224788348\n')
        self.assertEqual(expected, fp.getvalue())

    def testEmptyCopy(self):
        """
        The emptyCopy method must return a new, empty database.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        db.addSubject(subject)
        newDb = db.emptyCopy()
        self.assertIs(None, db.params.compare(newDb.params))
        self.assertEqual(0, newDb.subjectCount)
        self.assertEqual(0, newDb.addSubject(AARead('id1', 'AAA')))


class TestBackend(TestCase):
    """
    Tests for the light.database.Backend class.
    """
    def testParametersAreStored(self):
        """
        The backend must call its super class so its parameters are stored.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        self.assertIs(params, be.params)

    def testInitialBackendIsEmpty(self):
        """
        The index must be empty if no reads have been added.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        self.assertEqual({}, be.d)

    def testAddSubjectReturnsCorrectResult(self):
        """
        If one subject is added, addSubject must return the index ('0') of the
        added subject, the correct number of covered residues, and the correct
        hash count.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        subject = AARead('id', 'FRRRFRRRF')
        coveredResidues, hashCount = be.addSubject(subject, '0')
        self.assertEqual(9, coveredResidues)
        self.assertEqual(0, hashCount)

    def testAddSameSubjectIndex(self):
        """
        If the same subject index is passed more tha once to addSubject, it
        must raise ValueError.
        """
        params = Parameters([], [])
        be = Backend(params)
        subject = AARead('id', 'FRRRFRRRF')
        be.addSubject(subject, '0')
        error = "^Subject index '0' has already been used\.$"
        self.assertRaisesRegex(ValueError, error, be.addSubject, subject, '0')

    def testAddSameSubjectIncreasesBackendSize(self):
        """
        If an identical subject is added multiple times, the backend size
        increases, because the backend does not detect duplicates, only the
        Database frontend does.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRF'), '0')
        self.assertEqual(1, be.subjectCount)
        be.addSubject(AARead('id', 'FRRRFRRRF'), '1')
        self.assertEqual(2, be.subjectCount)

    def testOneReadOneLandmark(self):
        """
        If one subject is added but it only has one landmark, nothing is added
        to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRF'), '0')
        self.assertEqual({}, be.d)

    def testOneReadTwoLandmarks(self):
        """
        If one subject is added and it has two landmarks, one key is added
        to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(
            AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [0]},
            },
            be.d)

    def testTwoReadsTwoLandmarksLimitZeroPairsPerLandmark(self):
        """
        If two identical reads are added, both with two landmarks, no keys
        will be added to the dictionary if limitPerLandmark is zero.
        """
        params = Parameters([AlphaHelix], [], limitPerLandmark=0)
        be = Backend(params)
        be.addSubject(
            AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        self.assertEqual({}, be.d)

    def testTwoReadsTwoLandmarksDifferentOffsets(self):
        """
        If two subjects are added, both with two landmarks separated by the
        same distance, only one key is added to the backend and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the backend since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(
            AARead('id1', 'AFRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [1], '1': [0]},
            },
            be.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one subject is added and it has one landmark and one peak, one pair
        is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

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
        params = Parameters([AlphaHelix], [Peaks], distanceBase=distanceBase)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        distance10 = str(scale(10, distanceBase))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkOnePeakNoTrigFinders(self):
        """
        If one subject is added and it has one landmark and one peak, but no
        trig finders are in use, nothing is added to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        self.assertEqual({}, be.d)

    def testOneReadOneLandmarkTwoPeaks(self):
        """
        If one subject is added and it has one landmark and two peaks, two
        pairs are added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksLimitOnePairPerLandmark(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        limit of one pair per landmarks is imposed, only one key is added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], limitPerLandmark=1)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksSevereMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe maximum distance is imposed, no keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=1)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        self.assertEqual({}, be.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        maximum distance is imposed that makes one of the peaks too far
        away, only one key is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksLargeMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, and a
        maximum distance is imposed that is greater than the distance to the
        peaks, two keys are added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=15)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksPermissiveMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        permissive minimum distance is imposed, all keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=1)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
                'A2:P:' + distance10: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but an
        intermediate minimum distance is imposed, only the key for the pair
        that exceeds the minimum distance is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=11)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [0]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksSevereMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe minimum distance is imposed, no keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=100)
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        self.assertEqual({}, be.d)

    def testMultipleSubjectOffsets(self):
        """
        If one subject is added and it has one landmark and one peak separated
        by 10 bases and then, later in the subject, the same pair with the
        same separation, one key must be added to the backend and it
        should have two subject offsets.  Note that minDistance and
        maxDistance are used to discard the matches some longer and shorter
        distance pairs that only have one subject offset (i.e., that only
        appear in the subject once).
        """
        seq = 'FRRRFRRRFASA'
        params = Parameters([AlphaHelix], [Peaks], minDistance=5,
                            maxDistance=10)
        be = Backend(params)
        be.addSubject(AARead('id', seq + seq), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [0, 12]},
            },
            be.d)

    def testTwoReadsTwoLandmarksSameOffsets(self):
        """
        If two identical reads are added, both with two landmarks at the same
        offsets, only one key is added to the backend and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the backend since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(
            AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [0], '1': [0]},
            },
            be.d)

    def testFindNoMatch(self):
        """
        A query against an empty backend must produce no results.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRR')
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashInOneLocation(self):
        """
        One matching subject with one matching hash (that occurs in one
        location) must be found correctly.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, 0.0, Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual(
            {
                '0': [
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                        'queryLandmarkOffsets': [1],
                        'queryTrigPointOffsets': [11],
                        'subjectLandmarkOffsets': [1],
                        'trigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                ],
            },
            matches)
        self.assertEqual(1, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashInTwoLocations(self):
        """
        One matching subject with one matching hash (that occurs in two
        locations) must be found correctly.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASAVVVVVVASAVVVASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual(
            {
                '0': [
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
            }, matches)
        self.assertEqual(14, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindWithNonMatchingHashes(self):
        """
        Non-matching hashes must be found correctly when storeFullAnalysis is
        passed to find() as True.
        """
        subject = AARead('subject', 'F')
        query = AARead('query', 'AFRRRFRRRFASAASAVV')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, True)

        self.assertEqual({}, matches)
        self.assertEqual(2, hashCount)
        self.assertEqual(
            {
                'A2:P:10': {
                    'landmark':
                        Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 1, 9, 2),
                    'landmarkOffsets': [1],
                    'trigPoint': TrigPoint(Peaks.NAME, Peaks.SYMBOL, 11),
                    'trigPointOffsets': [11],
                },
                'A2:P:13': {
                    'landmark':
                        Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 1, 9, 2),
                    'landmarkOffsets': [1],
                    'trigPoint': TrigPoint(Peaks.NAME, Peaks.SYMBOL, 14),
                    'trigPointOffsets': [14],
                }
            },
            nonMatchingHashes)

    def testFindNoneMatchingTooSmallDistance(self):
        """
        No matches should be found if the max distance is too small.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=1)
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindNoneMatchingNoTrigPoint(self):
        """
        No matches should be found if there is only one landmark and there are
        no trig point finders.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindTwoMatchingInSameSubject(self):
        """
        Two matching hashes in the subject must be found correctly.
        """
        sequence = 'FRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)

        self.assertEqual(
            {
                '0': [
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
            matches)
        self.assertEqual(2, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testChecksumEmptyBackend(self):
        """
        The backend checksum must be the same as the checksum for its
        parameters when no subjects have been added to the backend.
        """
        params = Parameters([], [])
        be = Backend(params)
        self.assertEqual(params.checksum, be.checksum)

    def testChecksumAfterSubjectAdded(self):
        """
        The backend checksum must be as expected when a subject has been
        added to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        be.addSubject(subject, '0')

        checksum = Checksum(params.checksum).update([
            'id',
            sequence,
            '0',  # Hash count.
        ]).checksum
        self.assertEqual(checksum, be.checksum)

    def testScan(self):
        """
        The scan method must return a scanned subject.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        be.addSubject(subject, '0')
        scannedSubject = be.scan(subject)
        self.assertIsInstance(scannedSubject, ScannedRead)

    def testGetScannedPairs(self):
        """
        The getSequencePairs method must return pairs of
        (landmark, trigPoints).
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        be = Backend(params)
        be.addSubject(subject, '0')
        scannedSubject = be.scan(subject)
        pairs = list(be.getScannedPairs(scannedSubject))
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
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        query = AARead('query', 'AAA')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashesWithOneLandmark(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets. The result should be
        empty if there is only one landmark in the read.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        query = AARead('query', 'FRRRFRRRF')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets.
        """
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        be = Backend(params)
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFASAASA')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
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

    def testRestoreInvalidJSON(self):
        """
        If a backend restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        error = '^Expected object or value$'
        self.assertRaisesRegex(ValueError, error, Backend.restore,
                               StringIO('xxx'))

    def testSaveLoadWithNonDefaultParameters(self):
        """
        When asked to save and then load a backend with non-default
        parameters, a backend with the correct parameters must result.
        """
        params = Parameters([], [], limitPerLandmark=16, maxDistance=17,
                            minDistance=18, distanceBase=19.0)
        be = Backend(params)
        fp = StringIO()
        be.save(fp)
        fp.seek(0)
        result = be.restore(fp)
        self.assertIs(None, params.compare(result.params))
        self.assertEqual(be._subjectIndices, result._subjectIndices)

    def testSaveLoadNonEmpty(self):
        """
        When asked to save and then load a non-empty backend, the correct
        backend must result.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be = Backend(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        fp = StringIO()
        be.save(fp)
        fp.seek(0)
        result = Backend.restore(fp)
        self.assertEqual(be.subjectCount, result.subjectCount)
        self.assertEqual(be.d, result.d)
        self.assertEqual(be.params.landmarkClasses,
                         result.params.landmarkClasses)
        self.assertEqual(be.params.limitPerLandmark,
                         result.params.limitPerLandmark)
        self.assertEqual(be.params.maxDistance, result.params.maxDistance)
        self.assertEqual(be.params.minDistance, result.params.minDistance)
        self.assertEqual(be.params.trigPointClasses,
                         result.params.trigPointClasses)
        self.assertEqual(be.checksum, result.checksum)
        self.assertEqual(be._subjectIndices, result._subjectIndices)

    def testChecksumAfterSaveLoad(self):
        """
        A backend that has a sequence added to it, which is then saved and
        then re-loaded, and then has a second sequence is added to it must have
        the same checksum as a backend that simply has the two sequences added
        to it without interruption.
        """
        seq1 = 'FRRRFRRRFASAASA'
        seq2 = 'MMMMMMMMMFRRRFR'
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be1 = Backend(params1)
        be1.addSubject(AARead('id1', seq1), '0')
        fp = StringIO()
        be1.save(fp)
        fp.seek(0)
        be1 = Backend.restore(fp)
        be1.addSubject(AARead('id2', seq2), '1')

        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be2 = Backend(params2)
        be2.addSubject(AARead('id1', seq1), '0')
        be2.addSubject(AARead('id2', seq2), '1')

        self.assertEqual(be1.checksum, be2.checksum)

    def testEmptyCopy(self):
        """
        The emptyCopy method must return a new, empty backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        be.addSubject(subject, '0')
        newBe = be.emptyCopy()
        self.assertIs(None, be.params.compare(newBe.params))
        self.assertEqual({}, newBe.d)
        self.assertEqual(0, newBe.subjectCount)

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        fp = StringIO()
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        be = Backend(params)
        be.addSubject(subject, '0')
        be.print_(fp)
        expected = (
            'Hash count: 3\n'
            'Checksum: 1144016651\n'
            'Subjects (with offsets) by hash:\n'
            '   A2:P:10\n'
            '    0 [0]\n'
            '   A2:T:4\n'
            '    0 [0]\n'
            '   A2:T:8\n'
            '    0 [0]\n'
            'Landmark symbol counts:\n'
            '  AlphaHelix (A2): 3\n'
            'Trig point symbol counts:\n'
            '  Peaks (P): 1\n'
            '  Troughs (T): 2\n')
        self.assertEqual(expected, fp.getvalue())


class TestSimpleConnector(TestCase):
    """
    Tests for the light.database.SimpleConnector class.
    """

    def testAddSubjectReturnsCorrectResult(self):
        """
        If one subject is added, addSubject must return the name of the
        backend, the index ('0') of the added subject, the correct number
        of covered residues, and the correct hash count.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend(params)
        sc = SimpleConnector(be)
        subject = AARead('id', 'FRRRFRRRF')
        coveredResidues, hashCount = sc.addSubject(subject, '0')
        self.assertEqual(9, coveredResidues)
        self.assertEqual(0, hashCount)

    def testFindNoMatch(self):
        """
        A query on a simple connector that has an empty backend must produce
        no results.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRR')
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend(params)
        sc = SimpleConnector(be)
        sc.addSubject(subject, '0')
        result = sc.find(
            query, Parameters.DEFAULT_SIGNIFICANCE_METHOD,
            Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)
        self.assertEqual(1, len(result))
        matches, hashCount, nonMatchingHashes = result[0]

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashInOneLocation(self):
        """
        One matching subject with one matching hash (that occurs in one
        location) must be found correctly.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend(params)
        sc = SimpleConnector(be)
        sc.addSubject(subject, '0')
        result = sc.find(
            query, 0.0, Parameters.DEFAULT_SCORE_METHOD,
            Parameters.DEFAULT_SIGNIFICANCE_FRACTION, False)
        self.assertEqual(1, len(result))
        matches, hashCount, nonMatchingHashes = result[0]

        self.assertEqual(
            {
                '0': [
                    {
                        'landmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                        'queryLandmarkOffsets': [1],
                        'queryTrigPointOffsets': [11],
                        'subjectLandmarkOffsets': [1],
                        'trigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                ],
            },
            matches)
        self.assertEqual(1, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        fp = StringIO()
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        be = Backend(params)
        sc = SimpleConnector(be)
        sc.addSubject(subject, '0')
        sc.print_(fp)
        expected = (
            'Backend \'localhost\':\n'
            'Hash count: 3\n'
            'Checksum: 1144016651\n'
            'Subjects (with offsets) by hash:\n'
            '   A2:P:10\n'
            '    0 [0]\n'
            '   A2:T:4\n'
            '    0 [0]\n'
            '   A2:T:8\n'
            '    0 [0]\n'
            'Landmark symbol counts:\n'
            '  AlphaHelix (A2): 3\n'
            'Trig point symbol counts:\n'
            '  Peaks (P): 1\n'
            '  Troughs (T): 2\n')
        self.assertEqual(expected, fp.getvalue())
