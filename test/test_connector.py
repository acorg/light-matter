from unittest import TestCase

from dark.reads import AARead

from light.exceptions import BackendException
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.parameters import Parameters
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
        params = Parameters([AlphaHelix], [])
        sc = SimpleConnector(params)
        subject = AARead('id', 'FRRRFRRRF')
        preExisting, subjectIndex, hashCount = sc.addSubject(subject)
        self.assertFalse(preExisting)
        self.assertEqual('0', subjectIndex)
        self.assertEqual(0, hashCount)

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        sc = SimpleConnector(params)
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        sc.addSubject(subject)
        expected = (
            'Backends:\n'
            '  Name: backend\n'
            '  Hash count: 3\n'
            '  Checksum: 20379718\n'
            '  Subjects (with offsets) by hash:\n'
            '    A2:P:10\n'
            '      0 [[0, 9, 10]]\n'
            '    A2:T:4\n'
            '      0 [[0, 9, 4]]\n'
            '    A2:T:8\n'
            '      0 [[0, 9, 8]]\n'
            '  Landmark symbol counts:\n'
            '    AlphaHelix (A2): 3\n'
            '  Trig point symbol counts:\n'
            '    Peaks (P): 1\n'
            '    Troughs (T): 2')
        self.assertEqual(expected, sc.print_(printHashes=True))


class UnusedConnectorTests:

    def testAddUnknownBackend(self):
        """
        If an unknown backend is added, a BackendException must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        error = "^Unknown backend 'name'\.$"
        self.assertRaisesRegex(BackendException, error, db.addBackend, 'name')

    def testReconnectSameNameBackend(self):
        """
        If a backend tries to connect but re-uses an existing backend name,
        a BackendException must be raised.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        name, checksum, params = db.addBackend()
        error = "^Backend %r is already connected\.$" % name
        self.assertRaisesRegex(BackendException, error, db.addBackend, name)

    def testAddNewBackend(self):
        """
        When a new backend is added, the returned parameters must be those of
        the database.
        """
        params1 = Parameters([AlphaHelix], [Peaks])
        db = Database(params1)
        name, checksum, params2 = db.addBackend()
        self.assertIs(params1, params2)

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

        params = Parameters([AlphaHelix], [])
        db = Database(params, NoBackendConnector)
        query = AARead('id', 'AAA')
        error = "^No backends available\.$"
        self.assertRaisesRegex(BackendException, error, db.find, query)

    def testFindWithOneUnreconnectedBackend(self):
        """
        If a database has one unreconnected backend, calling find() must raise
        BackendException with the expected error message.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.disconnectedBackends['dummy'] = None
        query = AARead('id', 'AAA')
        error = "^Backend 'dummy' has not reconnected\.$"
        self.assertRaisesRegex(BackendException, error, db.find, query)

    def testFindWithTwoUnreconnectedBackends(self):
        """
        If a database has two unreconnected backends, calling find() must raise
        BackendException with the expected error message.
        """
        params = Parameters([AlphaHelix], [])
        db = Database(params)
        db.disconnectedBackends['dummy1'] = None
        db.disconnectedBackends['dummy2'] = None
        query = AARead('id', 'AAA')
        error = "^2 backends \(dummy1, dummy2\) have not reconnected\.$"
        self.assertRaisesRegex(BackendException, error, db.find, query)

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

        params = Parameters([AlphaHelix], [])
        db = Database(params, NoBackendConnector)
        error = "^Database has no backends\.$"
        self.assertRaisesRegex(BackendException, error, db.addSubject, None)

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

        params = Parameters([AlphaHelix], [])
        db = Database(params, TwoBackendConnector)
        self.assertEqual(['backend-0', 'backend-1'],
                         sorted(db.backends.keys()))
        self.assertEqual(Checksum().update('backend-0').checksum,
                         db.backends['backend-0'].checksum)
        self.assertEqual(Checksum().update('backend-1').checksum,
                         db.backends['backend-1'].checksum)
