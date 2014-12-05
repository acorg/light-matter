from light.features import TrigPoint


class AminoAcids(object):
    """
    A class for computing trig points based on the position of a specific
    amino acid.
    """
    NAME = 'AminoAcids'
    SYMBOL = 'M'

    def find(self, read, aa='C'):
        """
        A function that checks if and where a specific amino acid occurs in a
        sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @param aa: The single letter code of an amino acid. The default is 'C',
            as cysteines do not get substituted very often.
        """
        aa = aa or 'C'

        for i, base in enumerate(read.sequence):
            if base == aa:
                yield TrigPoint(self.NAME, self.SYMBOL, i)
