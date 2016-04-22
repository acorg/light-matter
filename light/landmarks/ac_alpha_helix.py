"""
Aho Corasick based alpha helix finder.

The Aho Corasick implementation is taken from
https://github.com/WojciechMula/pyahocorasick
"""

from os.path import dirname, join

import ahocorasick

import light
from light.features import Landmark
from light.finder import Finder


def _loadDatabase():
    """
    Read alpha helix prefix strings and add them to an Aho Corasick matcher.

    @return: An C{ahocorasick.Automaton} instance.
    """
    filename = join(dirname(light.__file__), '..', 'data',
                    'aho-corasick-alpha-helix-prefixes')
    # We only need to store the lengths of the prefixes, not the prefixes
    # themselves.
    ac = ahocorasick.Automaton(ahocorasick.STORE_LENGTH)
    add = ac.add_word
    with open(filename) as fp:
        for line in fp:
            add(line[:-1])
    ac.make_automaton()
    return ac


# Read the alpha helix database (only once) into a singleton for use in the
# ACAlphaHelix.find method.
_AC = _loadDatabase()


class AC_AlphaHelix(Finder):
    """
    A class for finding alpha helices using Aho Corasick to store prefixes.
    """
    NAME = 'AC AlphaHelix'
    SYMBOL = 'ACAH'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) alpha
        helices.

        No attempt is made to disallow overlapping alpha helices.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for end, length in _AC.iter(read.sequence):
            yield Landmark(self.NAME, self.SYMBOL, end - length + 1, length)