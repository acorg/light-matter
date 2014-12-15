from dark.aa import PROPERTY_DETAILS
from light.features import TrigPoint


class PolarityPeaks(object):
    """
    A class for computing statistics based on amino acid property peaks.
    The statistic sums up the values of all properties as it walks along a
    window on the amino acid sequence. In each window a peak is found, which
    has the highest aggregates property value.
    """
    NAME = 'PolarityPeak'
    SYMBOL = 'O'

    def convertAAToProperties(self, sequence, prop='polarity'):
        """
        Takes an amino acid sequence, converts it to a sequence of properties.

        @param sequence: an amino acid sequence.
        @param properties: a list of properties that should be included.
        """
        result = []
        previousProp = 0
        for aa in sequence:
            if aa in PROPERTY_DETAILS:
                aaProperty = PROPERTY_DETAILS[aa][prop] + previousProp
                result.append(aaProperty)
                previousProp += PROPERTY_DETAILS[aa][prop]
        return result

    def find(self, read, windowSize=50, prop='polarity'):
        """
        A function that checks if and where a peak helix in a sequence
        occurs.

        @param read: An instance of C{dark.reads.AARead}.
        @param properties: a list of properties that should be included.
        """
        for i in range(len(read.sequence) - windowSize):
            seq = read.sequence[i:i + 50]
            propertySequence = self.convertAAToProperties(seq, prop=prop)
            peakIndex = propertySequence.index(max(propertySequence)) + i
            yield TrigPoint(self.NAME, self.SYMBOL, peakIndex)
