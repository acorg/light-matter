from unittest import TestCase

import ahocorasick

from dark.reads import AARead

from light.features import Landmark
import light.landmarks.ac_alpha_helix
from light.landmarks.ac_alpha_helix import AC_AlphaHelix
from light.parameters import DatabaseParameters


def setAlphaHelices(helices):
    """
    Make an Aho Corasick matcher for the given helices and monkey patch
    light.landmarks.ac_alpha_helix to use it.
    Also set the acAlphaHelixFilename database parameter to 'xxx', to make sure
    that the right Aho Corasick matcher is used.

    This function is used by tests that want to check against a specific
    set of helices instead of the full set.

    @param helices: An interable of C{str} helix sequences.

    @return: A C{light.landmarks.ac_alpha_helix} instance with its
        acAlphaHelixFilename set to 'xxx'.
    """
    ac = ahocorasick.Automaton(ahocorasick.STORE_LENGTH)

    if ahocorasick.unicode:
        add = ac.add_word
    else:
        # The Aho Corasick module has been compiled without unicode
        # support, to reduce its memory use. Arrange to give it bytes.
        def add(s):
            """
            Add a string to an Aho Corasick automaton as bytes.

            @param s: A C{str} to add.
            """
            ac.add_word(s.encode('utf-8'))

    list(map(add, helices))
    ac.make_automaton()
    light.landmarks.ac_alpha_helix._AC = ac
    dbParams = DatabaseParameters(acAlphaHelixFilename='xxx')
    return AC_AlphaHelix(dbParams)


class TestACAlphaHelix(TestCase):
    """
    Tests for the light.landmarks.ac_alpha_helix.AC_AlphaHelix class.
    """
    def setUp(self):
        """
        Keep track of the original Aho Corasick matcher and stored filename so
        we can restore them after each test. This allows tests to install their
        own matcher and filename without breaking things for tests that want to
        use the original.
        """
        self.originalAC = light.landmarks.ac_alpha_helix._AC
        self.originalFile = light.landmarks.ac_alpha_helix._STORED_AC_FILENAME

    def tearDown(self):
        """
        Restore the original Aho Corasick state machine.
        """
        light.landmarks.ac_alpha_helix._AC = self.originalAC
        light.landmarks.ac_alpha_helix._STORED_AC_FILENAME = self.originalFile

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
        (in data/ac-alpha-helix-substrings-20-0.9).
        """
        read = AARead('id', 'RCELARTLKR VAWRN')
        finder = AC_AlphaHelix()
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelix', 'ACAH', 0, 4),
                Landmark('AC AlphaHelix', 'ACAH', 0, 5),
                Landmark('AC AlphaHelix', 'ACAH', 1, 4),
                Landmark('AC AlphaHelix', 'ACAH', 1, 5),
                Landmark('AC AlphaHelix', 'ACAH', 2, 6),
                Landmark('AC AlphaHelix', 'ACAH', 3, 6),
                Landmark('AC AlphaHelix', 'ACAH', 4, 5),
                Landmark('AC AlphaHelix', 'ACAH', 4, 6),
                Landmark('AC AlphaHelix', 'ACAH', 11, 4),
                Landmark('AC AlphaHelix', 'ACAH', 11, 5),
            ],
            sorted(result))
