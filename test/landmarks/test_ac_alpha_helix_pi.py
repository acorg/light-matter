from unittest import TestCase

import ahocorasick

from dark.reads import AARead

from light.features import Landmark
import light.landmarks.ac_alpha_helix_pi
from light.landmarks.ac_alpha_helix_pi import AC_AlphaHelix_pi
from light.parameters import DatabaseParameters


def setAlphaHelices_pi(helices):
    """
    Make an Aho Corasick matcher for the given helices and monkey patch
    light.landmarks.ac_alpha_helix_pi to use it.
    Also set the acAlphaHelixPiFilename database parameter to 'xxx', to make
    sure that the right Aho Corasick matcher is used.

    This function is used by tests that want to check against a specific
    set of helices instead of the full set.

    @param helices: An interable of C{str} helix sequences.

    @return: A C{light.landmarks.ac_alpha_helix_pi} instance with its
        acAlphaHelixPiFilename set to 'xxx'.
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
    light.landmarks.ac_alpha_helix_pi._AC = ac
    dbParams = DatabaseParameters(acAlphaHelixPiFilename='xxx')
    return AC_AlphaHelix_pi(dbParams)


class TestACAlphaHelix_pi(TestCase):
    """
    Tests for the light.landmarks.ac_alpha_helix_pi.AC_AlphaHelix_pi class.
    """
    def setUp(self):
        """
        Keep track of the original Aho Corasick matcher and stored filename so
        we can restore them after each test. This allows tests to install their
        own matcher and filename without breaking things for tests that want to
        use the original.
        """
        self.originalAC = light.landmarks.ac_alpha_helix_pi._AC
        self.originalFile = \
            light.landmarks.ac_alpha_helix_pi._STORED_AC_FILENAME

    def tearDown(self):
        """
        Restore the original Aho Corasick state machine.
        """
        light.landmarks.ac_alpha_helix_pi._AC = self.originalAC
        light.landmarks.ac_alpha_helix_pi._STORED_AC_FILENAME = \
            self.originalFile

    def testFindNothing(self):
        """
        The find method must return an empty generator when no helix is
        present.
        """
        finder = setAlphaHelices_pi(['XXX', 'YYY'])
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFullMatch(self):
        """
        The find method must return the full read sequence when it fully
        matches an alpha helix.
        """
        finder = setAlphaHelices_pi(['FFFF'])
        read = AARead('id', 'FFFF')
        result = list(finder.find(read))
        self.assertEqual([Landmark('AC AlphaHelixPi', 'ACAHP', 0, 4)], result)

    def testFindContiguousMatches(self):
        """
        The find method must find matches that are contiguous.
        """
        finder = setAlphaHelices_pi(['RRR', 'FFF'])
        read = AARead('id', 'FFFRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelixPi', 'ACAHP', 0, 3),
                Landmark('AC AlphaHelixPi', 'ACAHP', 3, 3),
            ],
            sorted(result))

    def testFindSeparatedMatches(self):
        """
        The find method must find matches that are separated.
        """
        finder = setAlphaHelices_pi(['RRRRR', 'FFF'])
        read = AARead('id', 'FFFMMRRRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelixPi', 'ACAHP', 0, 3),
                Landmark('AC AlphaHelixPi', 'ACAHP', 5, 5),
            ],
            sorted(result))

    def testFindPartiallyOverlappingMatches(self):
        """
        The find method must return overlapping helices.
        """
        finder = setAlphaHelices_pi(['FFFFR', 'FRMMM'])
        read = AARead('id', 'FFFFRMMM')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelixPi', 'ACAHP', 0, 5),
                Landmark('AC AlphaHelixPi', 'ACAHP', 3, 5),
            ],
            sorted(result))

    def testFindCompletelyOverlappingMatches(self):
        """
        The find method must return all helices, including those that overlap.
        """
        finder = setAlphaHelices_pi(['FF', 'FFF'])
        read = AARead('id', 'FFF')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelixPi', 'ACAHP', 0, 2),
                Landmark('AC AlphaHelixPi', 'ACAHP', 0, 3),
                Landmark('AC AlphaHelixPi', 'ACAHP', 1, 2),
            ],
            sorted(result))

    def testFindUsingBuiltInAlphaHelices(self):
        """
        The find method must be able to find helices in the default alpha
        helix substring file loaded by light/landmarks/ac_alpha_helix_pi.py
        (in data/ac-alpha-helix-pi-substrings-1-0.5).
        """
        read = AARead('id', 'RCELARTLKRLREGIG')
        finder = AC_AlphaHelix_pi()
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC AlphaHelixPi', 'ACAHP', 10, 6),
            ],
            sorted(result))
