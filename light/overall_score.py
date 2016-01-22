from light.string import MultilineString


class BestBinScore(object):
    """
    Assigns the score of the best bin to be the score of the overall match.

    @param histogram: A C{light.histogram} instance.
    @param significantBins: A C{list} of C{dict}'s where each dict contains
        information about the score, bin and index of a significant bin. This
        list is already sorted by score.
    """
    def __init__(self, histogram, significantBins):
        try:
            score = significantBins[0]['score']
        except IndexError:
            score = None
        self._score = score
        self._analysis = {
            'score': score,
            'scoreClass': self.__class__,
        }

    def calculateScore(self):
        """
        Calculates the overall score for all histogram bins.

        @return: a C{float} score for the best significant bin and a C{dict}
            with information about the score.
        """
        return self._score, self._analysis

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
            'Overall score method: %s' % analysis['scoreClass'].__name__,
            'Overall score: %s' % analysis['score'],
        ])

        return str(result)

ALL_OVERALL_SCORE_CLASSES = [BestBinScore]
