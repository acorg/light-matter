import re

from light.statistics import Statistic

ORTHOBUNYA_3 = 'TCATCACATGA'
ORTHOBUNYA_5 = 'TCGTGTGATGA'
HANTA_3 = 'ATCATCATCTG'
HANTA_5 = 'ATGATGAT'
NAIRO_3 = 'AGAGTTTCT'
NAIRO_5 = 'AGAAACTCT'
PHLEBO_3 = 'TGTGTTTC'
PHLEBO_5 = 'GAAACACA'
TOSPO_3 = 'TCTCGTTAG'
TOSPO_5 = 'CTAACGAGA'
ORTHOBUNYA_3_RC = 'TCATGTGATGA'
ORTHOBUNYA_5_RC = 'TCATCACACGA'
HANTA_5_RC = 'ATCATCAT'

ALL_BUNYAVIRIDAE = [ORTHOBUNYA_3, ORTHOBUNYA_5, ORTHOBUNYA_3_RC,
                    ORTHOBUNYA_5_RC, HANTA_3, HANTA_5, HANTA_5_RC, NAIRO_3,
                    NAIRO_5, PHLEBO_3, PHLEBO_5, TOSPO_3, TOSPO_5]

BUNYAVIRIDAES = re.compile('.*' + '|'.join(ALL_BUNYAVIRIDAE))


class Bunyaviridae(Statistic):
    """
    A class that looks for sequences specific to bunyaviridae.
    See Fields Virology p. 1246 for details.
    """
    NAME = 'Bunyaviridae'
    SUPPORTED_TYPES = 'dna'
    MIN_LENGTH = min(ALL_BUNYAVIRIDAE)

    def findBunyaviridae(self, sequence):
        """
        Checks for Bunyaviridae sequences in the sequence.

        @param sequence: a nucleotide sequence.
        """
        return BUNYAVIRIDAES.search(sequence)

    def _evaluate(self, read):
        """
        Evaluates whether the read contains any of the bunyavirus sequences.

        @param read: the sequence that should be checked.
        @return: C{bool} True if any of the bunyavirus sequences match.
        """
        bunyavirus = self.findBunyaviridae(read.sequence)
        if bunyavirus:
            return True
