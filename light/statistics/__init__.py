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
    return count


def find(worksOn='all'):
    """
    A function to find all subclasses that should be executed. To be used to
    find the statistics to be run.

    TODO: Needs to be updated with more specific filtering methods.
          Expand globals!

    @param worksOn: The type of class that should be found. Can eigher be
    'all', 'dna', 'protein'.
    """
    from light.statistics.all_bases import HasAllBases
    ALL_STATISTICS = [HasAllBases]

    if worksOn == 'all':
        return ALL_STATISTICS
    elif worksOn == 'dna':
        dnaStatistics = []
        for statistic in ALL_STATISTICS:
            if statistic.SUPPORTED_TYPES == 'dna':
                dnaStatistics.append(statistic)
        return dnaStatistics
    elif worksOn == 'protein':
        proteinStatistics = []
        for statistic in ALL_STATISTICS:
            if statistic.SUPPORTED_TYPES == 'protein':
                dnaStatistics.append(statistic)
        return proteinStatistics


# Default exports for 'from light.statistics import *'
__all__ = ['Statistic', 'runStatistic', 'find']
