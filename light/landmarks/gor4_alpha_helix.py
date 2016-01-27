from light.distance import scale
from light.features import Landmark
from light.finder import Finder


class GOR4AlphaHelix(Finder):
    """
    Use the GOR IV algorithm to do probabilistic identification of amino acid
    sequences that are (hopefully!) likely to correspond to Alpha helices.

    See the src/gor4 directory in the dark-matter repo at
    https://github.com/acorg/dark-matter for information about our
    GOR IV implementation.
    """
    NAME = 'GOR4AlphaHelix'
    SYMBOL = 'GA'

    def find(self, read):
        """
        Find possible alpha helices in a sequence, using GOR IV. Alpha helices
        are indicated by 'H' characters in the GOR IV prediction string.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        predictions = read.gor4()['predictions']
        length = 0
        for offset, prediction in enumerate(predictions):
            if prediction == 'H':
                if length:
                    # We're already in a string of H's. Keep counting.
                    length += 1
                else:
                    start = offset
                    length = 1
            else:
                if length:
                    # We were in a string of H's, but it has just ended.
                    yield Landmark(self.NAME, self.SYMBOL, start, length,
                                   scale(length, self._featureLengthBase))
                    length = 0

        if length:
            # We reached the end of the string still in an alpha helix.
            yield Landmark(self.NAME, self.SYMBOL, start, length,
                           scale(length, self._featureLengthBase))
