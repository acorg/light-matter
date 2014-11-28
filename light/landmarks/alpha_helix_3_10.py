import re

from dark.aa import PROPERTIES, HYDROPHOBIC
from light.features import Landmark

ALPHA_HELIX_3_10 = re.compile('OIIO(?:IIO)+')


class AlphaHelix_3_10(object):
    """
    A class for computing statistics based on 3-10 alpha helices.  Based
    around the assumption that a 3-10 alpha helix is composed of at least two
    repeats of one hydrophobic and then 2 hydrophilic amino acids.

    """
    NAME = 'AlphaHelix_3_10'
    SYMBOL = 'B'

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
        for match in ALPHA_HELIX_3_10.finditer(hhProperties):
            start = match.start()
            end = match.end()
            repeatCount = (end - start - 1) / 3
            yield Landmark(self.SYMBOL, start, end - start, repeatCount)
