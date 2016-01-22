from operator import itemgetter

from light.string import MultilineString


class OverallScore(object):
    """
    Calculates the overall score that takes into account all significant
    histogram bins.

    @param histogram: A C{light.histogram} instance.
    @param significantBins: A C{list} of C{dict}'s where each dict contains
        information about the score, bin and index of a significant bin.
    """
    def __init__(self, histogram, significantBins):
        self._histogram = histogram
        self._significantBins = significantBins
        self.scoreGetter = itemgetter('score')

    def calculateScore(self):
        """
        Calculates the overall score for all histogram bins.

        @return: a C{float} score for all significant bins and a C{dict} with
            information about the score.
        """
        self._significantBins.sort(key=self.scoreGetter, reverse=True)
        try:
            score = self._significantBins[0]['score']
        except IndexError:
            score = 0.0

        analysis = {'score': score}

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin=''):
        """
        Convert an analysis to a nicely formatted string.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} human-readable version of the last analysis.
        """
        result = MultilineString(margin=margin)

        result.extend([
            'Overall score',
            'Score: %(score).4f' % analysis,
        ])

        return str(result)
