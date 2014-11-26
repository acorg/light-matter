from unittest import TestCase

from light import database


class FakeCursor(object):
    def __init__(self, results):
        self._results = results
        self._index = -1

    def execute(self, p):
        pass

    def fetchone(self):
        self._index += 1
        return self._results[self._index]

    def close(self):
        pass


class FakeDbConnection(object):
    def __init__(self, results):
        self._results = results
        self.open = True

    def cursor(self):
        return FakeCursor(self._results)

    def close(self):
        self.open = False


class TestDatabase(TestCase):
    """
    Tests for the database functions.
    """
    def testCompareVirusPresent(self):
        """
        If the alpha helix motif is present in the database, it must be
        returned correctly.
        """
        db = FakeDbConnection([('Luckyvirus'), ('[0, 14]')])
        result = database.compare('0|14', db=db)
        self.assertEqual('Luckyvirus', result)

    def testCompareVirusNotPresent(self):
        """
        If the alpha helix motif is not present in the database, it must be
        returned correctly.
        """
        db = FakeDbConnection([None])
        result = database.compare([2, 20], db=db)
        self.assertEqual(None, result)
