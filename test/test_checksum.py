from unittest import TestCase

from light.checksum import Checksum


class TestChecksum(TestCase):
    """
    Tests for the light.database.Checksum class.
    """
    def testUnused(self):
        """
        A checksum that has not been given any text must have the default
        0x0 value.
        """
        self.assertEqual(0x0, Checksum().value)

    def testSetInitialValue(self):
        """
        A checksum that is initialized with a specific value but which has not
        been given any text must be set to the passed initial value.
        """
        self.assertEqual(1, Checksum(1).value)

    def testEquality(self):
        """
        Two checksums that have been given the same text must compare as equal.
        """
        c1 = Checksum().update('text')
        c2 = Checksum().update('text')
        self.assertEqual(c1, c2)

    def testInequality(self):
        """
        Two checksums that have been given different text must not (normally)
        compare as equal.
        """
        c1 = Checksum().update('text1')
        c2 = Checksum().update('text2')
        self.assertNotEqual(c1, c2)

    def testSet(self):
        """
        It must be possible to set the checksum's value.
        """
        c = Checksum()
        c.value = 10
        self.assertEqual(10, c.value)
