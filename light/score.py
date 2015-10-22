from warnings import warn
from operator import itemgetter
from copy import copy


class MinHashesScore(object):
    """
    Calculates the score for a bin based on the minimum hashes of either
    subject or query.

    @param histogram: A C{light.histogram} instance.
    @param minHashCount: The C{int} minimum number of hashes in the query
        and the subject.
    """
    def __init__(self, histogram, minHashCount):
        self._histogram = histogram
        self.minHashCount = minHashCount

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{float} of the score of that bin.
        """
        binCount = len(self._histogram[binIndex])
        try:
            score = float(binCount) / self.minHashCount
        except ZeroDivisionError:
            score = 0.0

        # It is possible that a bin has more deltas in it than
        # the number of hashes in the query / subject. This can
        # occur when the distanceBase used to scale distances
        # results in multiple different distances being put
        # into the same bin. In this case, the distance binning
        # is perhaps too aggessive. For now, issue a warning
        # and set the score to 1.0
        if score > 1.0:
            score = 1.0
            warn('Bin contains %d deltas for a query/subject pair '
                 'with a minimum hash count of only %d. The '
                 'database distanceBase is causing different '
                 'distances to be scaled into the same bin, which '
                 'might not be a good thing.' %
                 (binCount, self.minHashCount), RuntimeWarning)

        return score


def histogramBinFeatures(bin_, queryOrSubject):
    """
    Extract all features (of type C{queryOrSubject} (a C{str}, either
    'query' or 'subject') from a bin and return them as a set, along with
    the min and max offsets of the features.

    @param bin_: A C{light.histogram.Histogram} bin.
    @param queryOrSubject: A C{str}, to indicate which features to extract,
        either 'query' or 'subject'.
    @raise KeyError: If C{queryOrSubject} is not 'query' or 'subject'.
    @return: A 2-tuple, containing 1) a C{set} of all features (landmarks and
        trig points), 2) a C{set} of the offsets of all features in the bin
        (this includes the start and end of all landmarks, plus the offset
        of all trig points).
    """
    allOffsets = set()
    allFeatures = set()
    # There is no error checking that queryOrSubject is 'query' or
    # 'subject' as the following item getter will raise a KeyError if it
    # cannot access the dict key in the bin element.
    offsetsListGetter = itemgetter(queryOrSubject + 'Offsets')
    for hash_ in bin_:
        offsetsList = offsetsListGetter(hash_)
        landmark = hash_.landmark
        trigPoint = hash_.trigPoint
        for offsets in offsetsList:
            # Copy the landmark, set its offset, add it to our results.
            thisLandmark = copy(landmark)
            thisLandmark.offset = offsets[0]
            allFeatures.add(thisLandmark)
            allOffsets.add(thisLandmark.offset)
            allOffsets.add(thisLandmark.offset + thisLandmark.length)
            # Copy the trig point, set its offset, add it to our results.
            thisTrigPoint = copy(trigPoint)
            thisTrigPoint.offset = offsets[1]
            allFeatures.add(thisTrigPoint)
            allOffsets.add(thisTrigPoint.offset)

    return allFeatures, allOffsets


class FeatureMatchingScore:
    """
    Calculates the score for histogram bins based on the features present
    within the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    """
    def __init__(self, histogram, query, subject):
        self._histogram = histogram
        self._query = query
        self._subject = subject

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{float} of the score of that bin.
        """
        return 0.0


ALL_SCORE_CLASSES = (MinHashesScore, FeatureMatchingScore)
