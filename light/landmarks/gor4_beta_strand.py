from light.distance import scale
from light.features import Landmark, Finder


class GOR4BetaStrand(Finder):
    """
    Use the GOR IV algorithm to do probabilistic identification of amino acid
    sequences that are (hopefully!) likely to correspond to Beta strands.

    See the src/gor4 directory in the dark-matter repo at
    https://github.com/acorg/dark-matter for information about our
    GOR IV implementation.
    """
    NAME = 'GOR4BetaStrand'
    SYMBOL = 'GB'

    def find(self, read):
        """
        Find possible beta strands in a sequence, using GOR IV. Beta strands
        are indicated by 'E' characters in the GOR IV prediction string.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        predictions = read.gor4()['predictions']
        count = 0
        for offset, prediction in enumerate(predictions):
            if prediction == 'E':
                if count:
                    # We're already in a string of E's. Keep counting.
                    count += 1
                else:
                    start = offset
                    count = 1
            else:
                if count:
                    # We were in a string of E's, but it has just ended.
                    yield Landmark(self.NAME, self.SYMBOL, start, count,
                                   scale(count, self._featureLengthBase))
                    count = 0

        if count:
            # We reached the end of the string still in a beta strand.
            yield Landmark(self.NAME, self.SYMBOL, start, count,
                           scale(count, self._featureLengthBase))
