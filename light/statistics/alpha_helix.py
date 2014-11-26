from light.statistics import Statistic

PROPERTIES = {'F': 'O', 'Y': 'O', 'W': 'O', 'H': 'O', 'K': 'O', 'T': 'O',
              'C': 'O', 'G': 'O', 'A': 'O', 'V': 'O', 'I': 'O', 'L': 'O',
              'M': 'O', 'P': 'I', 'S': 'I', 'N': 'I', 'D': 'I', 'Q': 'I',
              'E': 'I', 'R': 'I', 'X': 'O'}

ALPHAHELIX = 'OIIIOIIIO'


class AlphaHelix(Statistic):
    """
    A class for computing statistics based on alpha helices.
    Based around the assumption that an alpha helix is composed of three times
    one hydrophobic and 3 hydrophilic amino acids.
    """
    NAME = 'AlphaHelix'
    SUPPORTED_TYPES = 'aa'
    MIN_LENGTH = 9

    def convertAAToHydrophobicHydrophilic(self, sequence):
        """
        Takes an amino acid sequence, converts it to a sequence of hydrophobic
        / hydrophilic.

        @param sequence: an amino acid sequence.
        """
        return ''.join(PROPERTIES[aa] for aa in sequence)

    def find(self, hhProperties):
        """
        A function that checks if and where an alpha helix in a sequence
        occurs.
        """
        positions = []
        for index in range(len(hhProperties) - 8):
            toEvaluate = hhProperties[index:index + 9]
            if toEvaluate == ALPHAHELIX:
                positions.append(index)
        return positions

    def calculateDistance(self, positions):
        """
        A function that takes a list of positions and calculates a list of
        distances between those positions. These are stored in a database to
        look up.

        @param positions: a C{List} of positions, as returned by find.
        """
        if len(positions) == 1:
            return False
        else:
            distances = []
            for i in range(len(positions) - 1):
                distance = positions[i + 1] - positions[i]
                distances.append(distance)
            return distances

    def _evaluate(self, read, distances=False):
        """
        Test if a sequence has an alpha helix.

        @param read: the sequence that should be checked.
        """
        readHHProperties = self.convertAAToHydrophobicHydrophilic(
            read.sequence)
        result = self.find(readHHProperties)
        if distances:
            result = self.calculateDistance(result)
            return result
        return len(result) > 0
