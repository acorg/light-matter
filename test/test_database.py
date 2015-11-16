from unittest import TestCase
from io import StringIO
try:
    from ujson import loads
except ImportError:
    from json import loads

from dark.reads import AARead

from light.features import Landmark, TrigPoint
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.checksum import Checksum
from light.parameters import Parameters, FindParameters
from light.database import Database
from light.subject import Subject
from light.backend import Backend


class TestDatabase(TestCase):
    """
    Tests for the light.database.Database class.
    """
    def testSignificanceFractionDefault(self):
        """
        The significanceFraction default value must be as expected.
        """
        self.assertEqual(0.25, FindParameters.DEFAULT_SIGNIFICANCE_FRACTION)

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
        self.assertEqual(0, db.subjectCount())
        self.assertEqual(0, db.totalResidues())
        self.assertEqual(0, db.totalCoveredResidues())

    def testInitialDatabaseHasNoSubjectInfo(self):
        """
        The database must not have any stored subject information if no
        subjects have been added.
        """
        params = Parameters([], [])
        db = Database(params)
        self.assertEqual([], list(db.getSubjects()))

    def testAddSubjectReturnsIndexAndPreExisting(self):
        """
        If one subject is added, addSubject must return the index ('0') of the
        added subject and a Boolean to indicate whether the subject was already
        in the database.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        preExisting, subjectIndex, hashCount = db.addSubject(
            AARead('id', 'FRRRFRRRF'))
        self.assertFalse(preExisting)
        self.assertEqual('0', subjectIndex)
        self.assertEqual(0, hashCount)

    def testAddSameSubjectReturnsSameIndex(self):
        """
        If an identical subject is added multiple times, the same subject
        index must be returned.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        self.assertEqual(db.addSubject(AARead('id', 'FRRRFRRRF'))[1],
                         db.addSubject(AARead('id', 'FRRRFRRRF'))[1])

    def testAddSameSubjectReturnsCorrectPreExisting(self):
        """
        If an identical subject is added multiple times, the expected
        pre-existing values must be returned.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        self.assertEqual([False, True],
                         [db.addSubject(AARead('id', 'FRRRFRRRF'))[0],
                          db.addSubject(AARead('id', 'FRRRFRRRF'))[0]])

    def testAddSameSubjectLeavesDatabaseSizeTheSame(self):
        """
        If an identical subject is added multiple times, the database size
        does not increase.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount())
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount())

    def testGetSubjectByIndexKeyError(self):
        """
        If an unknown subject index is passed to getSubjectByIndex, a KeyError
        must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        error = "^'xxx'$"
        self.assertRaisesRegex(KeyError, error, db.getSubjectByIndex, 'xxx')

    def testGetIndexBySubjectKeyError(self):
        """
        If a non-existent subject is passed to getIndexBySubject, a KeyError
        must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        error = "^'id'$"
        self.assertRaisesRegex(KeyError, error, db.getIndexBySubject,
                               Subject('id', 'FF', 5))

    def testGetSubjectByIndex(self):
        """
        If a subject is added, getSubjectByIndex must be able to return it
        given its string index.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        subject = AARead('id', 'FRRRFRRRF')
        _, index, _ = db.addSubject(subject)
        self.assertEqual(Subject('id', 'FRRRFRRRF', 0),
                         db.getSubjectByIndex(index))

    def testGetSubjectBySubject(self):
        """
        If a subject is added, getIndexBySubject must be able to return it
        given an identical subject to look up.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        _, index, _ = db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(index,
                         db.getIndexBySubject(Subject('id', 'FRRRFRRRF', 0)))

    def testGetSubjectHashCount(self):
        """
        If a subject is added, getSubjectByIndex must return a Subject
        instance that has the correct hash count.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        subject = AARead('id', 'FRRRFRRRFAFRRRFRRRF')
        _, index, _ = db.addSubject(subject)
        self.assertEqual(1, db.getSubjectByIndex(index).hashCount)

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
        self.assertEqual(Subject('id', 'FRRRFRRRF', 0), subject)
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
        self.assertEqual(Subject('id', 'FRRRFRRRFAFRRRFRRRF', 1), subject)
        self.assertEqual(1, subject.hashCount)

    def testOneReadOneLandmarkStatistics(self):
        """
        If one subject is added the database statistics must be correct.
        """
        params = Parameters([], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRF'))
        self.assertEqual(1, db.subjectCount())
        self.assertEqual(9, db.totalResidues())
        self.assertEqual(0, db.totalCoveredResidues())

    def testOneReadTwoLandmarksStatistics(self):
        """
        If one subject is added, the database statistics must be correct.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(1, db.subjectCount())
        self.assertEqual(36, db.totalResidues())
        self.assertEqual(22, db.totalCoveredResidues())

    def testTwoReadsTwoLandmarksStatistics(self):
        """
        If two identical reads are added, the database statistics must be
        correct.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.addSubject(AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        db.addSubject(AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'))
        self.assertEqual(2, db.subjectCount())
        self.assertEqual(72, db.totalResidues())
        self.assertEqual(44, db.totalCoveredResidues())

    def testSaveContentHasFourParts(self):
        """
        When a simple database saves, its content must include parts for the
        database parameters, the database state, the backend parameters and
        the backend state.
        """
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        dbParams = Parameters.restore(fp)
        loads(fp.readline()[:-1])
        backendParams = Parameters.restore(fp)
        loads(fp.readline()[:-1])
        self.assertIs(None, dbParams.compare(backendParams))

    def testSaveContentIncludesExpectedKeysAndValues(self):
        """
        When the database saves, its JSON content must include the expected
        keys and values.
        """
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        Parameters.restore(fp)
        state = loads(fp.readline()[:-1])

        # Keys
        self.assertEqual(['_connectorClassName'], list(state.keys()))

        # Values
        self.assertEqual('SimpleConnector', state['_connectorClassName'])

    def testSaveRestoreEmpty(self):
        """
        When asked to save and then restore an empty database, the correct
        database must result.
        """
        params = Parameters([AlphaHelix], [Peaks])
        db = Database(params)
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.restore(fp)
        self.assertEqual(0, result.subjectCount())
        self.assertEqual(0, result.totalCoveredResidues())
        self.assertEqual(0, result.totalResidues())
        self.assertIs(None, params.compare(result.params))

    def testRestoreInvalidJSON(self):
        """
        If a database restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        params = Parameters([], [Peaks])
        db = Database(params)
        error = '^Expected object or value$'
        self.assertRaisesRegex(ValueError, error, db.restore, StringIO('xxx'))

    def testSaveRestoreNonEmpty(self):
        """
        When asked to save and then restore a non-empty database, the correct
        database must result.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        db = Database(params)
        db.addSubject(AARead('id', 'FRRRFRRRFASAASA'))
        fp = StringIO()
        db.save(fp)
        fp.seek(0)
        result = db.restore(fp)
        self.assertIs(None, params.compare(result.params))
        self.assertEqual(db.subjectCount(), result.subjectCount())
        self.assertEqual(db.totalCoveredResidues(),
                         result.totalCoveredResidues())
        self.assertEqual(db.totalResidues(), result.totalResidues())
        self.assertEqual(db.checksum(), result.checksum())

    def testChecksumAfterSaveRestore(self):
        """
        A database that has a sequence added to it, which is then saved and
        restored, and then has a second sequence is added to it must have the
        same checksum as a database that simply has the two sequences added to
        it without the intervening save/restore.
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

        self.assertEqual(db1.checksum(), db2.checksum())

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

    def testFindMatchAfterSaveRestore(self):
        """
        A matching subject found before a save/restore must also be found
        following a database save/restore.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASAVVVVVVASAVVVASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db1 = Database(params)
        db1.addSubject(subject)
        result = db1.find(query)
        expected = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 10),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 13),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 14),
                }
            ]
        }
        self.assertEqual(expected, result.matches)
        fp = StringIO()
        db1.save(fp)
        fp.seek(0)
        db2 = Database.restore(fp)
        result = db2.find(query)
        self.assertEqual(expected, result.matches)

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
                '0': [
                    {
                        'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryTrigPoint': TrigPoint('Peaks', 'P', 10),
                        'subjectLandmark': Landmark(
                            'AlphaHelix', 'A', 1, 9, 2),
                        'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                    {
                        'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                        'queryTrigPoint': TrigPoint('Peaks', 'P', 13),
                        'subjectLandmark': Landmark(
                            'AlphaHelix', 'A', 1, 9, 2),
                        'subjectTrigPoint': TrigPoint('Peaks', 'P', 14),
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
        findParams = FindParameters(significanceFraction=0.0)
        result = db.find(query, findParams)
        self.assertEqual(
            {
                '0': [
                    {
                        'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                        'queryTrigPoint': TrigPoint('Peaks', 'P', 11),
                        'subjectLandmark': Landmark(
                            'AlphaHelix', 'A', 1, 9, 2),
                        'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                ],
            },
            result.matches)

    def testFindOneMatchingSignificantWithSubjectIndicesIncludingIt(self):
        """
        One matching and significant subject must be found, including when a
        non-empty subjectIndices is passed which includes the found index (and
        other non-matched subject indices)
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        db = Database(params)
        db.addSubject(subject)
        findParams = FindParameters(significanceFraction=0.0)
        result = db.find(query, findParams, subjectIndices={'0', 'x', 'y'})
        self.assertEqual(
            {
                '0': [
                    {
                        'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                        'queryTrigPoint': TrigPoint('Peaks', 'P', 11),
                        'subjectLandmark': Landmark(
                            'AlphaHelix', 'A', 1, 9, 2),
                        'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                    },
                ],
            },
            result.matches)

    def testFindOneMatchingButSubjectExcluded(self):
        """
        Despite one matching and significant subject, no result should be
        returned if a subjectIndices argument that excludes it is passed to
        find.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        db = Database(params)
        db.addSubject(subject)
        findParams = FindParameters(significanceFraction=0.0)
        result = db.find(query, findParams, subjectIndices=set())
        self.assertEqual({}, result.matches)

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
            {'0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 10),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 10),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 13),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 13),
                }
            ]
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
        findParams = FindParameters(significanceFraction=0.0)
        result = db.find(query, findParams)
        score1 = result.analysis['0']['bestScore']

        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params)
        db.addSubject(query)
        result = db.find(subject, findParams)
        score2 = result.analysis['0']['bestScore']

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
        _, index, _ = db.addSubject(subject)
        hashCount1 = db.getSubjectByIndex(index).hashCount
        findParams = FindParameters(significanceFraction=0.0)
        result = db.find(query, findParams)
        score1 = result.analysis['0']['bestScore']

        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks])
        db = Database(params2)
        _, index, _ = db.addSubject(query)
        hashCount2 = db.getSubjectByIndex(index).hashCount
        result = db.find(subject, findParams)
        score2 = result.analysis['0']['bestScore']

        self.assertNotEqual(hashCount1, hashCount2)
        self.assertEqual(score1, score2)
        self.assertNotEqual(1.0, score1)

    def testChecksumEmptyDatabase(self):
        """
        The database checksum must be the same as the checksum for its
        parameters plus the default backend name when no subjects have
        been added to the database.
        """
        params = Parameters([], [])
        expected = Checksum(params.checksum).update([
            Backend.DEFAULT_NAME,
        ])
        db = Database(params)
        self.assertEqual(expected.value, db.checksum())

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

        expected = Checksum(params.checksum).update([
            Backend.DEFAULT_NAME,
            'id',
            sequence,
        ])
        self.assertEqual(expected.value, db.checksum())

    def testSaveRestoreWithNonDefaultParameters(self):
        """
        When asked to save and then restore a database with non-default
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

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
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
            '  Distance base: 1.000000\n'
            'Connector class: SimpleConnector\n'
            'Subject count: 1\n'
            'Hash count: 3\n'
            'Total residues: 15\n'
            'Coverage: 73.33%\n'
            'Checksum: 539152060\n'
            'Connector:')
        self.assertEqual(expected, db.print_())

    def testPrintWithHashes(self):
        """
        The print_ function should produce the expected output when asked to
        print hash information.
        """
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
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
            '  Distance base: 1.000000\n'
            'Connector class: SimpleConnector\n'
            'Subject count: 1\n'
            'Hash count: 3\n'
            'Total residues: 15\n'
            'Coverage: 73.33%\n'
            'Checksum: 20379718\n'
            'Connector:\n'
            'Backends:\n'
            '  Name: backend\n'
            '  Hash count: 3\n'
            '  Checksum: 20379718\n'
            '  Subjects (with offsets) by hash:\n'
            '    A2:P:10\n'
            '      0 [[0, 9, 10, 1]]\n'
            '    A2:T:4\n'
            '      0 [[0, 9, 4, 1]]\n'
            '    A2:T:8\n'
            '      0 [[0, 9, 8, 1]]\n'
            '  Landmark symbol counts:\n'
            '    AlphaHelix (A2): 3\n'
            '  Trig point symbol counts:\n'
            '    Peaks (P): 1\n'
            '    Troughs (T): 2')
        self.assertEqual(expected, db.print_(printHashes=True))

    def testPrintNoHashes(self):
        """
        The print_ function should report the expected result if no hashes are
        found in the subject.
        """
        subject = AARead('subject', '')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        db = Database(params)
        db.addSubject(subject)
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
            '  Distance base: 1.000000\n'
            'Connector class: SimpleConnector\n'
            'Subject count: 1\n'
            'Hash count: 0\n'
            'Total residues: 0\n'
            'Coverage: 0.00%\n'
            'Checksum: 425336937\n'
            'Connector:')
        self.assertEqual(expected, db.print_())
