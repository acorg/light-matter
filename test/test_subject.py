from unittest import TestCase

from dark.reads import AARead

from light.subject import Subject


class TestSubject(TestCase):
    """
    Tests for the light.database.Subject class.
    """
    def testHashCountIsStored(self):
        """
        A Subject must store the hashcount it is passed.
        """
        self.assertEqual(6, Subject('id', 'AA', 6).hashCount)

    def testIsAARead(self):
        """
        A Subject is a subclass of AARead.
        """
        self.assertTrue(isinstance(Subject('id', 'AA', 0), AARead))

    def testAAReadInitCalled(self):
        """
        A Subject must call AARead.__init__ with the correct arguments. We can
        check this by comparing an AARead to a Subject (and vice versa).
        """
        self.assertEqual(AARead('id', 'A', '!'), Subject('id', 'A', 6, '!'))
        self.assertEqual(Subject('id', 'A', 6, '!'), AARead('id', 'A', '!'))
