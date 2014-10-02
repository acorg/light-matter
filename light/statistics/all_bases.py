from light.statistics import Statistic


class HasAllBases(Statistic):
    """
    Test that DNA sequences contain all four nucleotides.
    """
    NAME = 'hasAllBases'
    SUPPORTED_TYPES = 'dna'
    MIN_LENGTH = 4

    def _evaluate(self, sequence):
        seq = sequence.sequence
        return 'A' in seq and 'C' in seq and 'T' in seq and 'G' in seq
