from unittest import TestCase

import ahocorasick

from dark.reads import AARead

from light.features import Landmark
import light.landmarks.ac_extended_strand
from light.landmarks import AC_ExtendedStrand
from light.parameters import DatabaseParameters


def setExtendedStrands(strands):
    """
    Make an Aho Corasick matcher for the given strands and monkey patch
    light.landmarks.ac_extended_strand to use it.
    Also set the acExtendedStrandFilename database parameter to 'xxx', to make
    sure that the right Aho Corasick matcher is used.

    This function is used by tests that want to check against a specific
    set of strands instead of the full set.

    @param strands: An iterable of C{str} strand sequences.

    @return: A C{light.landmarks.ac_extended_strand} instance with its
        acExtendedStrandFilename set to 'xxx'.
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

    list(map(add, strands))
    ac.make_automaton()
    light.landmarks.ac_extended_strand._AC = ac
    dbParams = DatabaseParameters(acExtendedStrandFilename='xxx')
    return AC_ExtendedStrand(dbParams)


class TestACExtendedStrand(TestCase):
    """
    Tests for the light.landmarks.ac_extended_strand.AC_ExtendedStrand class.
    """
    def setUp(self):
        """
        Keep track of the original Aho Corasick matcher and stored filename so
        we can restore them after each test. This allows tests to install their
        own matcher and filename without breaking things for tests that want to
        use the original.
        """
        self.originalAC = light.landmarks.ac_extended_strand._AC
        self.originalFile = \
            light.landmarks.ac_extended_strand._STORED_AC_FILENAME

    def tearDown(self):
        """
        Restore the original Aho Corasick state machine.
        """
        light.landmarks.ac_extended_strand._AC = self.originalAC
        light.landmarks.ac_extended_strand._STORED_AC_FILENAME = \
            self.originalFile

    def testFindNothing(self):
        """
        The find method must return an empty generator when no strand is
        present.
        """
        finder = setExtendedStrands(['XXX', 'YYY'])
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFullMatch(self):
        """
        The find method must return the full read sequence when it fully
        matches an extended strand.
        """
        finder = setExtendedStrands(['FFFF'])
        read = AARead('id', 'FFFF')
        result = list(finder.find(read))
        self.assertEqual([Landmark('AC ExtendedStrand', 'ACES', 0, 4)], result)

    def testFindContiguousMatches(self):
        """
        The find method must find matches that are contiguous.
        """
        finder = setExtendedStrands(['RRR', 'FFF'])
        read = AARead('id', 'FFFRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC ExtendedStrand', 'ACES', 0, 3),
                Landmark('AC ExtendedStrand', 'ACES', 3, 3),
            ],
            sorted(result))

    def testFindSeparatedMatches(self):
        """
        The find method must find matches that are separated.
        """
        finder = setExtendedStrands(['RRRRR', 'FFF'])
        read = AARead('id', 'FFFMMRRRRR')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC ExtendedStrand', 'ACES', 0, 3),
                Landmark('AC ExtendedStrand', 'ACES', 5, 5),
            ],
            sorted(result))

    def testFindPartiallyOverlappingMatches(self):
        """
        The find method must return overlapping strands.
        """
        finder = setExtendedStrands(['FFFFR', 'FRMMM'])
        read = AARead('id', 'FFFFRMMM')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC ExtendedStrand', 'ACES', 0, 5),
                Landmark('AC ExtendedStrand', 'ACES', 3, 5),
            ],
            sorted(result))

    def testFindCompletelyOverlappingMatches(self):
        """
        The find method must return all strands, including those that overlap.
        """
        finder = setExtendedStrands(['FF', 'FFF'])
        read = AARead('id', 'FFF')
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC ExtendedStrand', 'ACES', 0, 2),
                Landmark('AC ExtendedStrand', 'ACES', 0, 3),
                Landmark('AC ExtendedStrand', 'ACES', 1, 2),
            ],
            sorted(result))

    def testFindUsingBuiltInExtendedStrands(self):
        """
        The find method must be able to find strands in the default extended
        strand substring file loaded by light/landmarks/ac_extended_strand.py
        (in data/ac-extended-strand-substrings-10-0.5).
        """
        read = AARead('id', 'RCELARTLKRFCCC')
        finder = AC_ExtendedStrand()
        result = list(finder.find(read))
        self.assertEqual(
            [
                Landmark('AC ExtendedStrand', 'ACES', 9, 4),
                Landmark('AC ExtendedStrand', 'ACES', 10, 4),
            ],
            sorted(result))
