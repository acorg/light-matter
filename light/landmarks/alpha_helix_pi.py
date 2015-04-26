import re

from dark.aa import PROPERTIES, HYDROPHOBIC
from light.features import Landmark, Finder

ALPHA_HELIX_PI = re.compile('OIIIIO(?:IIIIO)+')


class AlphaHelix_pi(Finder):
    """
    A class for computing statistics based on Pi alpha helices.  Based
    around the assumption that a pi alpha helix is composed of at least two
    repeats of one hydrophobic and then 4 hydrophilic amino acids.

    """
    NAME = 'AlphaHelix_pi'
    SYMBOL = 'C'

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
        for match in ALPHA_HELIX_PI.finditer(hhProperties):
            start = match.start()
            end = match.end()
            symbolDetail = (end - start - 1) // 5
            yield Landmark(self.NAME, self.SYMBOL, start, end - start,
                           symbolDetail)
