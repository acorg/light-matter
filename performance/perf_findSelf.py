from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks import BetaStrand
from light.trig import Peaks


class _FindSelfMixin(object):
    """
    Test if a database build from a sequence finds itself.
    """
    SEQUENCE = ''
    def testFindSelf(self):
        """
        Does a sequence match itself using different landmark and trig point
        finders.
        """



class TestFindSelf(TestCase):
    """
    Test that a database built from one sequence can match the same sequence.
    """

    def testFindIdenticalSequence(self):
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
