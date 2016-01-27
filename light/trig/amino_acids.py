from light.features import TrigPoint
from light.finder import Finder


class AminoAcids(Finder):
    """
    A class for computing trig points based on the position of a specific
    amino acid.
    """
    NAME = 'AminoAcids'
    SYMBOL = 'M'

    def find(self, read, aa=None):
        """
        A function that checks if and where a specific amino acid occurs in a
        sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @param aa: A list of the single letter codes of amino acids. The
            default is [W'].
        """
        aa = aa or ['W']

        for i, base in enumerate(read.sequence):
            if base in aa:
                yield TrigPoint(self.NAME, self.SYMBOL, i)
