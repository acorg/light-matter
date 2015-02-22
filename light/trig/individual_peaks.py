from dark.aa import PROPERTY_DETAILS as PD
from light.features import TrigPoint, Finder


class IndividualPeaks(Finder):
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
            middleAA = read.sequence[i]
            afterAA = read.sequence[i + 1]
            if beforeAA in PD and middleAA in PD and afterAA in PD:
                before = {'composition': PD[beforeAA]['composition'],
                          'iep': PD[beforeAA]['iep'],
                          'polarity': PD[beforeAA]['polarity']}

                middle = {'composition': PD[middleAA]['composition'],
                          'iep': PD[middleAA]['iep'],
                          'polarity': PD[middleAA]['polarity']}

                after = {'composition': PD[afterAA]['composition'],
                         'iep': PD[afterAA]['iep'],
                         'polarity': PD[afterAA]['polarity']}
                if (before['composition'] < middle['composition'] and
                        after['composition'] < middle['composition'] and
                        before['iep'] < middle['iep'] and
                        after['iep'] < middle['iep'] and
                        before['polarity'] < middle['polarity'] and
                        after['polarity'] < middle['polarity']):
                    yield TrigPoint(self.NAME, self.SYMBOL, i)
