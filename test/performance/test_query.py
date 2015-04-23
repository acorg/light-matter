from unittest import TestCase

from dark.reads import Reads, AARead

from light.database import Database
from light.performance.query import queryDatabase
from light.landmarks import AlphaHelix
from light.trig import Peaks


class TestQueryDatabase(TestCase):
    """
    Tests for the queryDatabase function.

    Note that we are not trying to test that database building or look-up
    works (that is all tested elsewhere), just that the queryDatabase
    function returns the expected structure with the expected result.
    """

    def testNoMatches(self):
        """
        No matches should be found if the database is empty.
        """
        database = Database([AlphaHelix], [Peaks], maxDistance=11)
        subjects = Reads()
        queries = Reads()
        queries.add(AARead('query', 'FRRRFRRRFASAASA'))
        result = queryDatabase(subjects, queries, database,
                               significanceMethod='hashFraction',
                               significanceFraction=0.1)
        self.assertEqual({}, result)

    def testFindOneMatchingSignificant(self):
        """
        One matching and significant subject must be found if the
        significanceFraction is sufficiently low.
        """
        database = Database([AlphaHelix], [Peaks], maxDistance=11)
        subjects = Reads()
        subjects.add(AARead('subject', 'AFRRRFRRRFASAASA'))
        queries = Reads()
        queries.add(AARead('query', 'FRRRFRRRFASAASA'))
        result = queryDatabase(subjects, queries, database,
                               significanceFraction=0.1)
        self.assertEqual(
            {
                'query': {
                    'subject': 1,
                },
            },
            result)
