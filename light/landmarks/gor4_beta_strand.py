from light.distance import scaleLog
from light.features import Landmark
from light.finder import Finder

from gor4 import GOR4


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
        gor4 = GOR4()
        predictions = gor4.predict(read.sequence)
        featureLengthBase = self._dbParams.featureLengthBase
        length = 0
        for offset, prediction in enumerate(predictions['predictions']):
            if prediction == 'E':
                if length:
                    # We're already in a string of E's. Keep counting.
                    length += 1
                else:
                    start = offset
                    length = 1
            else:
                if length:
                    # We were in a string of E's, but it has just ended.
                    yield Landmark(self.NAME, self.SYMBOL, start, length,
                                   scaleLog(length, featureLengthBase))
                    length = 0

        if length:
            # We reached the end of the string still in a beta strand.
            yield Landmark(self.NAME, self.SYMBOL, start, length,
                           scaleLog(length, featureLengthBase))
