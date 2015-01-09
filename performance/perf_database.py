import time
from unittest import TestCase

from dark.reads import AARead

from light.database import Database


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

    def testAdd10KSubjects(self):
        """
        How long does it take to add 10K subjects to a database that has no
        landmark or trig point finders?
        """
        read = AARead('id', 'MTMTSTTSNLPGILSQPSSELLTNWYAEQVVQGHIL')
        database = Database([], [])
        startTime = time.time()
        for _ in xrange(10000):
            database.addSubject(read)
        elapsed = time.time() - startTime
        self.details = elapsed
