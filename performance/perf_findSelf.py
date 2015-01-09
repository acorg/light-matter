from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks import BetaStrand
from light.trig import Peaks


class TestFindSelf(TestCase):
    """
    Test that a database built from one sequence can match the same sequence.
    """

    def testFindIdenticalSequenced(self):
        """
        Build a database with one subject and then search for that same
        sequence. Add details of the result to self.
        """
        seq = 'MTMTSTTSNLPGILSQPSSELLTNWYAEQVVQGHIL'
        read = AARead('id', seq)
        database = Database([BetaStrand], [Peaks])
        database.addSubject(read)
        result = database.find(read)
        significant = list(result.significant())
        if 0 in significant:
            self.details = {
                'result': True,
                'score': result.analysis[0]['score'],
            }
        else:
            self.details = {
                'result': False,
            }
