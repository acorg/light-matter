import six
from unittest import TestCase

from dark.reads import AARead

from light.exceptions import BackendException
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.parameters import DatabaseParameters
from light.database import Database
from light.connector import SimpleConnector
from light.checksum import Checksum


class TestSimpleConnector(TestCase):
    """
    Tests for the light.database.SimpleConnector class.
    """

    def testAddSubjectReturnsCorrectResult(self):
        """
        If one subject is added, addSubject must return the correct
        pre-existing status and subject index.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        sc = SimpleConnector(dbParams)
        subject = AARead('id', 'FRRRFRRRF')
        preExisting, subjectIndex, hashCount = sc.addSubject(subject)
        self.assertFalse(preExisting)
        self.assertEqual('0', subjectIndex)
        self.assertEqual(0, hashCount)

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix, BetaStrand],
                                      trigPoints=[Peaks, Troughs],
                                      limitPerLandmark=16, maxDistance=10,
                                      minDistance=0, distanceBase=1)
        sc = SimpleConnector(dbParams)
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        sc.addSubject(subject)
        expected = (
            'Backends:\n'
            '  Name: backend\n'
            '  Hash count: 3\n'
            '  Checksum: 337886368\n'
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
        self.assertEqual(expected, sc.print_(printHashes=True))


class UnusedConnectorTests(object):

    def testAddUnknownBackend(self):
        """
        If an unknown backend is added, a BackendException must be raised.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams)
        error = "^Unknown backend 'name'\.$"
        six.assertRaisesRegex(self, BackendException, error, db.addBackend,
                              'name')

    def testReconnectSameNameBackend(self):
        """
        If a backend tries to connect but re-uses an existing backend name,
        a BackendException must be raised.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams)
        name, checksum, dbParams = db.addBackend()
        error = "^Backend %r is already connected\.$" % name
        six.assertRaisesRegex(self, BackendException, error, db.addBackend,
                              name)

    def testAddNewBackend(self):
        """
        When a new backend is added, the returned parameters must be those of
        the database.
        """
        dbParams1 = DatabaseParameters(landmarks=[AlphaHelix],
                                       trigPoints=[Peaks])
        db = Database(dbParams1)
        name, checksum, dbParams2 = db.addBackend()
        self.assertIs(dbParams1, dbParams2)

    def testFindWithNoBackends(self):
        """
        If no backend has been added to a database, calling find() must raise
        BackendException.
        """
        class NoBackendConnector(object):
            """
            A connector that does not add a backend to its database.

            @param database: A C{Database} instance.
            """
            def __init__(self, database):
                pass

        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams, NoBackendConnector)
        query = AARead('id', 'AAA')
        error = "^No backends available\.$"
        six.assertRaisesRegex(self, BackendException, error, db.find, query)

    def testFindWithOneUnreconnectedBackend(self):
        """
        If a database has one unreconnected backend, calling find() must raise
        BackendException with the expected error message.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams)
        db.disconnectedBackends['dummy'] = None
        query = AARead('id', 'AAA')
        error = "^Backend 'dummy' has not reconnected\.$"
        six.assertRaisesRegex(self, BackendException, error, db.find, query)

    def testFindWithTwoUnreconnectedBackends(self):
        """
        If a database has two unreconnected backends, calling find() must raise
        BackendException with the expected error message.
        """
        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams)
        db.disconnectedBackends['dummy1'] = None
        db.disconnectedBackends['dummy2'] = None
        query = AARead('id', 'AAA')
        error = "^2 backends \(dummy1, dummy2\) have not reconnected\.$"
        six.assertRaisesRegex(self, BackendException, error, db.find, query)

    def testAddSubjectWithNoBackends(self):
        """
        If no backend has been added to a database, calling addSubject() must
        raise BackendException.
        """
        class NoBackendConnector(object):
            """
            A connector that does not add a backend to its database.

            @param database: A C{Database} instance.
            """
            def __init__(self, database):
                pass

        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams, NoBackendConnector)
        error = "^Database has no backends\.$"
        six.assertRaisesRegex(self, BackendException, error, db.addSubject,
                              None)

    def testTwoBackends(self):
        """
        If two backends are added to a database, they must have the expected
        names and checksums.
        """
        class TwoBackendConnector(object):
            """
            A connector that adds two backends to its database.

            @param database: A C{Database} instance.
            """
            def __init__(self, database):
                database.addBackend()
                database.addBackend()

        dbParams = DatabaseParameters(landmarks=[AlphaHelix])
        db = Database(dbParams, TwoBackendConnector)
        self.assertEqual(['backend-0', 'backend-1'],
                         sorted(db.backends.keys()))
        self.assertEqual(Checksum().update('backend-0').checksum,
                         db.backends['backend-0'].checksum)
        self.assertEqual(Checksum().update('backend-1').checksum,
                         db.backends['backend-1'].checksum)
