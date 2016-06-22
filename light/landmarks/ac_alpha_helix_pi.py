"""
Aho Corasick based pi alpha helix finder.

The Aho Corasick implementation is taken from
https://github.com/WojciechMula/pyahocorasick
"""

import ahocorasick

from light.features import Landmark
from light.finder import Finder


def _loadDatabase(filename):
    """
    Read pi alpha helix strings and add them to an Aho Corasick matcher.

    @param filename: A C{str} filename of a file with pi alpha helices.

    @return: An C{ahocorasick.Automaton} instance.
    """
    # We only need to store the lengths of the helices, not the helices
    # themselves.
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

    with open(filename) as fp:
        for line in fp:
            add(line[:-1])
    ac.make_automaton()
    return ac


# Read the pi alpha helix database (only once) into a singleton for use in the
# AC_AlphaHelix.find method.
_AC = None

# Store the ac alpha helix pi filename for checking.
# For testing, we want to prevent the if statement in find() from running.
# Therefore, we set _STORED_AC_FILENAME to 'xxx', so that
# self._dbParams.acAlphaHelixPiFilename can also be set to 'xxx'.
_STORED_AC_FILENAME = 'xxx'


class AC_AlphaHelix_pi(Finder):
    """
    A class for finding pi alpha helices using Aho Corasick to store helices.
    """
    NAME = 'AC AlphaHelixPi'
    SYMBOL = 'ACAHP'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) pi alpha
        helices.
        No attempt is made to disallow overlapping pi alpha helices.

        @param read: An instance of C{dark.reads.AARead}.

        @return: A generator that yields C{Landmark} instances.
        """
        global _AC, _STORED_AC_FILENAME
        if _AC is None or (
                _STORED_AC_FILENAME != self._dbParams.acAlphaHelixPiFilename):
            _AC = _loadDatabase(self._dbParams.acAlphaHelixPiFilename)
            _STORED_AC_FILENAME = self._dbParams.acAlphaHelixPiFilename

        if ahocorasick.unicode:
            sequence = read.sequence
        else:
            sequence = read.sequence.encode('utf-8')

        for end, length in _AC.iter(sequence):
            yield Landmark(self.NAME, self.SYMBOL, end - length + 1, length)
