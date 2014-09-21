class Statistic(object):

    def evaluate(self, sequence):
        """
        @param sequence: a dark.Read object
        """
        try:
            self.SUPPORTED_TYPES
        except AttributeError:
            raise Exception('Has no attribute SUPPORTED_TYPES')

        try:
            self.MIN_LENGTH
        except AttributeError:
            raise Exception('Has no attribute MIN_LENGTH')

        if sequence.type not in self.SUPPORTED_TYPES:
            return False

        if len(sequence) < self.MIN_LENGTH:
            return False

        return self._evaluate(sequence)


class HasAllBases(Statistic):
    """
    Test that DNA sequences contain all four nucleotides.
    """
    NAME = 'hasAllBases'
    SUPPORTED_TYPES = ['dna']
    MIN_LENGTH = 4

    def _evaluate(self, sequence):
        seq = sequence.sequence
        return 'A' in seq and 'C' in seq and 'T' in seq and 'G' in seq
