import time
from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.beta_strand import BetaStrand


class TestDatabase(TestCase):
    """
    Test database build performance.
    """

    def testCreation(self):
        """
        How long does it take to construct a database with just one subject if
        there are no landmark or trig point finders?
        """
        database = Database([], [])
        read = AARead('id', 'MTMTSTTSNLPGILSQPSSELLTNWYAEQVVQGHIL')
        database.addSubject(read)

    def testChecksumEmpty(self):
        """
        How long does it take to compute a checksum on an empty database?
        """
        database = Database([], [])
        startTime = time.time()
        database.checksum()
        elapsed = time.time() - startTime
        self.details = elapsed

    def testChecksum10K(self):
        """
        How long does it take to compute a checksum on a database that has 10K
        sequences added to it?
        """
        read = AARead('id', 'MTMTSTTSNLPGILSQPSSELLTNWYAEQVVQGHIL')
        database = Database([], [])
        for _ in xrange(10000):
            database.addSubject(read)
        startTime = time.time()
        database.checksum()
        elapsed = time.time() - startTime
        self.details = elapsed


class _TestSubjectAddMixin(object):
    """
    Test database subject adding performance.
    """

    def testAdd1000Subjects(self):
        """
        How long does it take to add 1000 subjects to a database.
        """
        read = AARead('id', 'MTMTSTTSNLPGILSQPSSELLTNWYAEQVVQGHIL')
        database = Database(self.LANDMARKS, self.TRIG_POINTS)
        startTime = time.time()
        for _ in xrange(1000):
            database.addSubject(read)
        elapsed = time.time() - startTime
        self.details = elapsed


class SubjectAddNoLandmarksNoTrigPoints(_TestSubjectAddMixin, TestCase):
    LANDMARKS = []
    TRIG_POINTS = []


class SubjectAddBetaStrand(_TestSubjectAddMixin, TestCase):
    LANDMARKS = [BetaStrand]
    TRIG_POINTS = []


class SubjectAddAlphaHelix(_TestSubjectAddMixin, TestCase):
    LANDMARKS = [AlphaHelix]
    TRIG_POINTS = []
