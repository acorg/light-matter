from dark.aa import PROPERTY_DETAILS
from light.features import TrigPoint


class Peaks(object):
    """
    A class for computing statistics based on amino acid property peaks.
    """
    NAME = 'Peaks'
    SYMBOL = 'P'

    def convertAAToProperties(self, sequence, properties=None):
        """
        Takes an amino acid sequence, converts it to a sequence of properties.

        @param sequence: an amino acid sequence.
        @param properties: a list of properties that should be included.
        """
        properties = properties or ['composition', 'iep', 'polarity']
        result = []

        for aa in sequence:
            aaProperties = sum(PROPERTY_DETAILS[aa][prop] for prop in
                               properties)
            result.append(aaProperties)
        return result

    def find(self, read, properties=None):
        """
        A function that checks if and where a peak helix in a sequence
        occurs.

        @param read: An instance of C{dark.reads.AARead}.
        @param properties: a list of properties that should be included.
        """
        properties = properties or ['composition', 'iep', 'polarity']
        propertySequence = self.convertAAToProperties(read.sequence,
                                                      properties)
        for i in range(1, len(read.sequence) - 1):
            before = propertySequence[i - 1]
            middle = propertySequence[i]
            after = propertySequence[i + 1]
            if before < middle and after < middle:
                yield TrigPoint(self.NAME, self.SYMBOL, i)