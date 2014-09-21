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


def runStatistic(statistic, fastaReads):
    """
    Function that runs the statistic on a set of reads.

    @param statistic: the name of a statistic class to be calculated.
    @param fastaReads: a dark.Reads object.
    """
    count = 0
    for read in fastaReads:
        result = statistic.evaluate(read)
        if result:
            count += 1
    return statistic.NAME, count


# Default exports for 'from light.statistics import *'
__all__ = ['Statistic', 'runStatistic']
