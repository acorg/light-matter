from light.features import Landmark


class AminoAcids(object):
    """
    A class for computing Landmarks based on the position of a specific
    amino acid.
    """
    NAME = 'AminoAcidsLm'
    SYMBOL = 'N'

    def find(self, read, aa=None):
        """
        A function that checks if and where a specific amino acid occurs in a
        sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @param aa: A list of the single letter codes of amino acids. The
            default is ['C'], as cysteines do not get substituted very
            often.
        """
        aa = aa or ['C']

        for i, base in enumerate(read.sequence):
            if base in aa:
                yield Landmark(self.NAME, self.SYMBOL, i, 1)
