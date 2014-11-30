from unittest import TestCase

from light.database import ScannedReadDatabase
from dark.reads import AARead


class TestScannedReadDatabase(TestCase):
    """
    Tests for the ScannedReadDatabase class.
    """
    def testInitialDictionary(self):
        """
        The dictionary must be instantiated correctly.
        """
        db = ScannedReadDatabase(['AlphaHelix'], ['Peaks'])
        self.assertEqual(['Peaks'],
                         db.d['params']['trigPointFinders'])
        self.assertEqual(['AlphaHelix'],
                         db.d['params']['landmarkFinders'])

    def testMakeSearchDictionary(self):
        """
        The search dictionary must be constructed correctly.
        """
        read = AARead('id', 'AFFFAAAFFFAAAA')
        db = ScannedReadDatabase(['AlphaHelix'],
                                 ['Peaks']).makeSearchDictionary(read)
        print db.d
        self.assertEqual('zzz', db.d)
