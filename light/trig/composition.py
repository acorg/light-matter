from dark.aa import PROPERTY_CLUSTERS
from light.features import TrigPoint, Finder


class Composition(Finder):
    """
    A class for computing trig points based on similarity in chemical
    composition.
    """
    NAME = 'Composition'
    SYMBOL = 'CP'

    def convertAAToCompositionClusters(self, sequence):
        """
        Takes an amino acid sequence, converts it to a sequence of numbers
        corresponding to the cluster an amino acid composition is in.

        @param sequence: an amino acid sequence.

        @return: A C{list} that for each amino acid contains the number of the
            cluster that amino acid is in.
        """
        return [PROPERTY_CLUSTERS[aa]['composition'] if aa in PROPERTY_CLUSTERS
                else 0 for aa in sequence]

    def find(self, read):
        """
        A function that checks if and where a change of composition in a
        sequence occurs.

        @param read: An instance of C{dark.reads.AARead}.
        """
        propertySequence = self.convertAAToCompositionClusters(read.sequence)
        for i in range(1, len(propertySequence) - 1):
            before = propertySequence[i - 1]
            middle = propertySequence[i]
            after = propertySequence[i + 1]
            if before == after and before != middle:
                yield TrigPoint(self.NAME, self.SYMBOL, i)
