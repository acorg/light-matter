from dark.aa import PROPERTY_DETAILS
from light.features import TrigPoint


class PolarityPeaks(object):
    """
    A class for computing statistics based on amino acid property peaks.
    The statistic sums up the values of all properties as it walks along a
    window on the amino acid sequence. In each window a peak is found, which
    has the highest aggregated property value.
    """
    NAME = 'PolarityPeak'
    SYMBOL = 'O'

    def sumProperties(self, sequence, prop='polarity'):
        """
        Takes an amino acid sequence, converts it to a sequence of properties.

        @param sequence: an amino acid sequence.
        @param prop: a C{str} name of the property which should be used to
            calculate the peak.
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
        A function that checks if and where a polarity peak in a sequence
        occurs.

        @param read: An instance of C{dark.reads.AARead}.
        @param windowSize: a C{int} of the size of the window in which the
            peak should be found.
        @param prop: a C{str} name of the property which should be used to
            calculate the peak.
        """
        sumPropertiesSequence = self.sumProperties(read.sequence, prop=prop)
        for i in range(len(sumPropertiesSequence) - windowSize):
            sumPropertiesWindow = sumPropertiesSequence[i:i + windowSize]
            peakIndex = sumPropertiesWindow.index(max(sumPropertiesWindow)) + i
            yield TrigPoint(self.NAME, self.SYMBOL, peakIndex)
