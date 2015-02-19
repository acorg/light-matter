from dark.aa import PROPERTY_DETAILS
from light.features import TrigPoint


class PolarityPeaks(object):
    """
    A class for computing statistics based on amino acid property peaks.
    The statistic computes a list where each element is the value of the
    property of the amino acid at that position plus the values of that
    property of all previous amino acids. It then moves a window along the
    computed list and returns the index of the highest value as the offset of a
    trig point, provided a trig point with the same absolute offset has not
    been yielded before.
    """
    NAME = 'PolarityPeak'
    SYMBOL = 'O'

    def sumProperties(self, sequence, prop='polarity'):
        """
        Takes an amino acid sequence, and returns a list where each element
        is the value of the property of the amino acid at that position, plus
        the values of the properties of all previous amino acids.

        @param sequence: a C{str} amino acid sequence.
        @param prop: a C{str} name of the property which should be used to
            calculate the peak.

        @return: a C{list}, where each element is the sum of the polarity value
            of the amino acid at that position and the polarity values of all
            previous amino acids.
        """
        result = []
        previousSumProperty = 0
        for aa in sequence:
            if aa in PROPERTY_DETAILS:
                currentSumProperty = (PROPERTY_DETAILS[aa][prop] +
                                      previousSumProperty)
                result.append(currentSumProperty)
                previousSumProperty += PROPERTY_DETAILS[aa][prop]
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
        storedOffset = None
        sumPropertiesSequence = self.sumProperties(read.sequence, prop=prop)
        for i in range(len(sumPropertiesSequence) - windowSize):
            sumPropertiesWindow = sumPropertiesSequence[i:i + windowSize]
            peakIndex = sumPropertiesWindow.index(max(sumPropertiesWindow)) + i
            # make it so that each polarity peak only returns one trig point,
            # even if it occurs in more than one window.
            if peakIndex - 1 != storedOffset and peakIndex != storedOffset:
                storedOffset = peakIndex
                yield TrigPoint(self.NAME, self.SYMBOL, peakIndex)
            else:
                storedOffset = peakIndex
