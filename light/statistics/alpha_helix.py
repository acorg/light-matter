from light.statistics import Statistic
from light.utils import convertAAToAAProperties, HYDROPHOBIC, HYDROPHILIC


ALPHAHELIX = 'OIIIOIIIO'


class AlphaHelix(Statistic):
    """
    A class for computing statistics based on alpha helices.
    Based around the assumption that an alpha helix is composed of three times
    one hydrophobic and 3 hydriphilic amino acids.
    """
    NAME = 'AlphaHelix'
    SUPPORTED_TYPES = 'protein'
    MIN_LENGTH = 9

    def convertAAToHydrophobicHydrophilic(sequence):
        """
        Takes an amino acid sequence, converts it to a sequence of hydrophobic
        / hydrophilic.

        @param sequence: an amino acid sequence.
        """
        hhProperties = convertAAToAAProperties(sequence,
                                               [HYDROPHOBIC, HYDROPHILIC])
        return str(hhProperties)

    def find(hhProperties):
        """
        A function that checks if and where an alpha helix in a sequence
        occurs.

        """
        positions = []
        for index in range(len(hhProperties) - 9):
            toEvaluate = hhProperties[index:index + 9]
            if toEvaluate == ALPHAHELIX:
                positions.append(index)
        return positions

    def _evaluate(read):
        """
        Test if a sequence has an alpha helix.

        @param read: the sequence that should be checked.
        """
        readHHProperties = AlphaHelix.convertAAToHydrophobicHydrophilic(read)
        print readHHProperties
        result = AlphaHelix.find(readHHProperties)
        if len(result) > 0:
            return True
