import re

from dark.aa import PROPERTIES, HYDROPHOBIC
from light.features import Landmark

ALPHAHELIX = re.compile('OIIIO(?:IIIO)+')


class AlphaHelix(object):
    """
    A class for computing statistics based on alpha helices.  Based
    around the assumption that an alpha helix is composed of three times
    one hydrophobic and 3 hydrophilic amino acids.

    """
    NAME = 'AlphaHelix'
    SYMBOL = 'A'

    def convertAAToHydrophobicHydrophilic(self, sequence):
        """
        Takes an amino acid sequence, converts it to a sequence of hydrophobic
        / hydrophilic.

        @param sequence: an amino acid sequence.
        """
        result = []
        append = result.append
        get = PROPERTIES.get

        for residue in sequence:
            # 'O' comes from hydrophObic, 'I' from hydrophIlic
            append('O' if (get(residue, 0x0) & HYDROPHOBIC) else 'I')
        return ''.join(result)

    def find(self, read):
        """
        A function that checks if and where an alpha helix in a sequence
        occurs.

        @param read: An instance of C{dark.reads.AARead}.
        """
        hhProperties = self.convertAAToHydrophobicHydrophilic(read.sequence)
        for match in ALPHAHELIX.finditer(hhProperties):
            start = match.start()
            end = match.end()
            repeatCount = (end - start - 1) >> 2
            yield Landmark(self.SYMBOL, start, end - start, repeatCount)
