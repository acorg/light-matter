from light.features import Landmark


class GOR4Coil(object):
    """
    Use the GOR IV algorithm to do probabilistic identification of amino acid
    sequences that are (hopefully!) likely to correspond to coils.

    See the src/gor4 directory in the dark-matter repo at
    https://github.com/acorg/dark-matter for information about our
    GOR IV implementation.
    """
    NAME = 'GOR4Coil'
    SYMBOL = 'GC'

    def find(self, read):
        """
        Find possible coils in a sequence, using GOR IV. Coils are indicated by
        'C' characters in the GOR IV prediction string.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        predictions = read.gor4()['predictions']
        count = 0
        for offset, prediction in enumerate(predictions):
            if prediction == 'C':
                if count:
                    # We're already in a string of C's. Keep counting.
                    count += 1
                else:
                    start = offset
                    count = 1
            else:
                if count:
                    # We were in a string of C's, but it has just ended.
                    yield Landmark(self.NAME, self.SYMBOL, start, count, count)
                    count = 0

        if count:
            # We reached the end of the string still in a coil.
            yield Landmark(self.NAME, self.SYMBOL, start, count, count)
