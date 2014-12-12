from dark.aa import PROPERTY_DETAILS
from light.features import TrigPoint


class IndividualPeaks(object):
    """
    A class for computing statistics based on individual amino acid property
    peaks.
    """
    NAME = 'IndividualPeaks'
    SYMBOL = 'I'

    def find(self, read, properties=None):
        """
        A function that checks if and where a peak helix in a sequence
        occurs.

        @param read: An instance of C{dark.reads.AARead}.
        @param properties: a list of properties that should be included.
        """
        properties = properties or ['composition', 'iep', 'polarity']

        for i in range(1, len(read.sequence) - 1):
            beforeAA = read.sequence[i - 1]
            before = {'composition': PROPERTY_DETAILS[beforeAA]['composition'],
                      'iep': PROPERTY_DETAILS[beforeAA]['iep'],
                      'polarity': PROPERTY_DETAILS[beforeAA]['polarity']}
            middleAA = read.sequence[i]
            middle = {'composition': PROPERTY_DETAILS[middleAA]['composition'],
                      'iep': PROPERTY_DETAILS[middleAA]['iep'],
                      'polarity': PROPERTY_DETAILS[middleAA]['polarity']}
            afterAA = read.sequence[i + 1]
            after = {'composition': PROPERTY_DETAILS[afterAA]['composition'],
                     'iep': PROPERTY_DETAILS[afterAA]['iep'],
                     'polarity': PROPERTY_DETAILS[afterAA]['polarity']}
            if (before['composition'] < middle['composition'] and
                    after['composition'] < middle['composition'] and
                    before['iep'] < middle['iep'] and
                    after['iep'] < middle['iep'] and
                    before['polarity'] < middle['polarity'] and
                    after['polarity'] < middle['polarity']):
                yield TrigPoint(self.NAME, self.SYMBOL, i)
