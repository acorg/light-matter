from unittest import TestCase

import ahocorasick

from dark.reads import AARead

from light.features import Landmark
import light.landmarks.ac_alpha_helix
from light.landmarks.ac_alpha_helix import AC_AlphaHelix
from light.parameters import DatabaseParameters

# Keep track of the original Aho Corasick matcher so we can restore it
# after each test. This allows tests to install their own matcher without
# breaking things for tests that want to use the original.
_ORIGINAL_AC = light.landmarks.ac_alpha_helix._AC


def setAlphaHelices(helices):
    """
    Make an Aho Corasick matcher for the given helices and monkey patch
    light.landmarks.ac_alpha_helix to use it.
    Also set the ahocorasickFilename database parameter to 'xxx', to make sure
    that the right Aho Corasick matcher is used.

    This function is used by tests that want to check against a specific
    set of helices instead of the full set.

    @param helices: An interable of C{str} helix sequences.

    @return: A C{light.landmarks.ac_alpha_helix} instance with its
        ahocorasickFilename set to 'xxx'.
    """
    ac = ahocorasick.Automaton(ahocorasick.STORE_LENGTH)
    list(map(ac.add_word, helices))
    ac.make_automaton()
    light.landmarks.ac_alpha_helix._AC = ac
    dbParams = DatabaseParameters(ahocorasickFilename='xxx')
    return AC_AlphaHelix(dbParams)


class TestACAlphaHelix(TestCase):
    """
    Tests for the light.landmarks.ac_alpha_helix.AC_AlphaHelix class.
    """

    def tearDown(self):
        """
        Restore the original Aho Corasick state machine.
        """
        light.landmarks.ac_alpha_helix._AC = _ORIGINAL_AC

    def testFindNothing(self):
        """
        The find method must return an empty generator when no helix is
        present.
        """
        finder = setAlphaHelices(['XXX', 'YYY'])
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFullMatch(self):
        """
        The find method must return the full read sequence when it fully
        matches an alpha helix.
        """
        finder = setAlphaHelices(['FFFF'])
        read = AARead('id', 'FFFF')
        result = list(finder.find(read))
        self.assertEqual([Landmark('AC AlphaHelix', 'ACAH', 0, 4)], result)

    def testFindContiguousMatches(self):
        """
        The find method must find matches that are contiguous.
        """
        finder = setAlphaHelices(['RRR', 'FFF'])
        read = AARead('id', 'FFFRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 3),
                Landmark('AC AlphaHelix', 'ACAH', 3, 3),
            ],
            sorted(result))

    def testFindSeparatedMatches(self):
        """
        The find method must find matches that are separated.
        """
        finder = setAlphaHelices(['RRRRR', 'FFF'])
        read = AARead('id', 'FFFMMRRRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 3),
                Landmark('AC AlphaHelix', 'ACAH', 5, 5),
            ],
            sorted(result))

    def testFindPartiallyOverlappingMatches(self):
        """
        The find method must return overlapping helices.
        """
        finder = setAlphaHelices(['FFFFR', 'FRMMM'])
        read = AARead('id', 'FFFFRMMM')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 5),
                Landmark('AC AlphaHelix', 'ACAH', 3, 5),
            ],
            sorted(result))

    def testFindCompletelyOverlappingMatches(self):
        """
        The find method must return all helices, including those that overlap.
        """
        finder = setAlphaHelices(['FF', 'FFF'])
        read = AARead('id', 'FFF')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 2),
                Landmark('AC AlphaHelix', 'ACAH', 0, 3),
                Landmark('AC AlphaHelix', 'ACAH', 1, 2),
            ],
            sorted(result))

    def testFindUsingBuiltInAlphaHelices(self):
        """
        The find method must be able to find helices in the default alpha
        helix prefix file loaded by light/landscapes/ac_alpha_helix.py
        (in data/aho-corasick-alpha-helix-prefixes).
        """
        read = AARead('id', 'RCELARTLKR VAWRN')
        finder = AC_AlphaHelix()
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 4),
                Landmark('AC AlphaHelix', 'ACAH', 0, 10),
                Landmark('AC AlphaHelix', 'ACAH', 11, 4),
                Landmark('AC AlphaHelix', 'ACAH', 11, 5),
            ],
            sorted(result))

        # Set _STORED_AC_FILENAME back to 'xxx' so that it doesn't interfere
        # with the following tests.
        light.landmarks.ac_alpha_helix._STORED_AC_FILENAME = 'xxx'
