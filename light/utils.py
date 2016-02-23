import sys
from itertools import chain

_None = object()

if sys.version_info < (3, 4):
    def maxWithDefault(a, default=_None):
        try:
            return max(a)
        except ValueError:
            if default is _None:
                raise
            else:
                return default

    def minWithDefault(a, default=_None):
        try:
            return min(a)
        except ValueError:
            if default is _None:
                raise
            else:
                return default
else:
    maxWithDefault = max
    minWithDefault = min


def stringSpans(s):
    """
    Find regions of a string that contain repeated characters.

    E.g., list(stringSpans('aabcccdd')) = [('a', 0, 2),
                                           ('b', 2, 3),
                                           ('c', 3, 6),
                                           ('d', 6, 8)]

    @param s: A C{str} to examine.
    @return: A generator that yields triples. Each triple contains 1) the
        character from C{s} that makes up the next region, 2) the offset
        of the beginning of the region, 3) the offset one beyond the end
        of this contiguous region.  I.e., the final two elements of each triple
        could be used Python-style to extract the substring from C{s} if that
        were necessary. The 2nd element can be subtracted from the 3rd to
        get the length of the region. See the above example.
    """
    previous = None
    for offset, symbol in enumerate(chain(s, [object()])):
        if previous is None:  # Start of string.
            previous = symbol
            startOffset = 0
        elif symbol != previous:
            yield (previous, startOffset, offset)
            previous = symbol
            startOffset = offset
