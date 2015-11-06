from copy import copy
from itertools import filterfalse
from warnings import warn


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
    Extract all features of type C{queryOrSubject} (a C{str}, either
    'query' or 'subject') from a bin and return them as a set, along with
    the offsets contained in all the features.

    @param bin_: A C{light.histogram.Histogram} bin.
    @param queryOrSubject: A C{str}, to indicate which features to extract,
        either 'query' or 'subject'.
    @raise KeyError: If C{queryOrSubject} is not 'query' or 'subject'.
    @return: A 2-tuple, containing 1) a C{set} of all features (landmarks and
        trig points) in the hashes in C{bin_}, 2) a C{set} of all offsets of
        all features (of type C{queryOrSubject}) in the bin.
    """
    allFeatures = set()
    allOffsets = set()
    # There is no error checking that queryOrSubject is 'query' or
    # 'subject' as the following will raise a KeyError if it cannot access
    # the dict key in the bin element.
    for hashInfo in bin_:
        for suffix in 'Landmark', 'TrigPoint':
            feature = hashInfo[queryOrSubject + suffix]
            allFeatures.add(feature)
            allOffsets.update(feature.coveredOffsets())

    return allFeatures, allOffsets


def featureInRange(feature, minOffset, maxOffset):
    """
    Does a feature fall within a min/max offset range?

    @param feature: A C{light.features._Feature} subclass (i.e., a
        landmark or a trig point).
    @param minOffset: The minimum allowed offset for the feature start.
    @param maxOffset: The maximum allowed offset for the feature end.
    @return: A C{bool} to indicate whether the feature falls completely
        within the allowed range.
    """
    return feature.offset >= minOffset and (
        (feature.offset + feature.length - 1) <= maxOffset)


def getHashFeatures(readHashes):
    """
    Extract all features from hashes returned by backend.getHashes().

    @param readHashes: The result from calling backend.getHashes().
    @return: A C{set} of all features (landmarks and trig points) that are in
        the hashes returned from calling backend.getHashes().
    """
    allFeatures = set()

    for hashInfo in readHashes.values():
        for offsetPair in hashInfo['offsets']:
            # Make sure all landmarks and trigPoints added to allFeatures have
            # the correct offsets. Note that a hash might occur in more than
            # one location on the sequence, which means that a landmark or
            # trigPoint will have to be added to allFeatures twice but with
            # different offsets.
            landmark = copy(hashInfo['landmark'])
            trigPoint = copy(hashInfo['trigPoint'])
            landmark.offset, trigPoint.offset = offsetPair

            allFeatures.update([landmark, trigPoint])

    return allFeatures


class FeatureMatchingScore:
    """
    Calculates the score for histogram bins based on the features present
    within the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param params: A C{Parameters} instance.
    @param findParams: An instance of C{light.parameters.FindParameters} or
        C{None} to use default find parameters.
    """
    def __init__(self, histogram, query, subject, params, findParams=None):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)
        from light.parameters import FindParameters
        self._findParams = findParams or FindParameters()
        from light.backend import Backend
        backend = Backend()
        backend.configure(params)
        scannedQuery = backend.scan(query)
        self._allQueryFeatures = set(scannedQuery.landmarks +
                                     scannedQuery.trigPoints)
        scannedSubject = backend.scan(subject)
        self._allSubjectFeatures = set(scannedSubject.landmarks +
                                       scannedSubject.trigPoints)

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{float} of the score of that bin.
        """
        queryFeatures, queryOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'query')
        subjectFeatures, subjectOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'subject')

        matchScore = self._findParams.featureMatchScore * (
            len(queryFeatures) + len(subjectFeatures))

        minQueryOffset = min(queryOffsets, default=0)
        maxQueryOffset = max(queryOffsets, default=self._queryLen)
        minSubjectOffset = min(subjectOffsets, default=0)
        maxSubjectOffset = max(subjectOffsets, default=self._subjectLen)

        # The mismatch score is applied to all features that are not
        # among those in the bin and which fall inside the max and min
        # offsets of the features in the bin.
        mismatchScore = self._findParams.featureMismatchScore * (
            len(list(filter(
                lambda f: featureInRange(f, minQueryOffset, maxQueryOffset),
                self._allQueryFeatures - queryFeatures))) +
            len(list(filter(
                lambda f: featureInRange(f, minSubjectOffset,
                                         maxSubjectOffset),
                self._allSubjectFeatures - subjectFeatures))))

        return matchScore + mismatchScore


class FeatureAAScore:
    """
    Calculates the score for histogram bins based on the count of amino acids
    present in the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param params: A C{Parameters} instance.
    """
    def __init__(self, histogram, query, subject, params):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)

        from light.backend import Backend
        backend = Backend()
        backend.configure(params)

        scannedQuery = backend.scan(query)
        allQueryHashes = backend.getHashes(scannedQuery)
        self._allQueryFeatures = getHashFeatures(allQueryHashes)
        allQueryFeaturesLength = len(scannedQuery.coveredIndices())

        scannedSubject = backend.scan(subject)
        allSubjectHashes = backend.getHashes(scannedSubject)
        self._allSubjectFeatures = getHashFeatures(allSubjectHashes)
        allSubjectFeaturesLength = len(scannedSubject.coveredIndices())

        # Work out whether the query or the subject have more covered indices.
        self._normaliseByQuery = True
        if allQueryFeaturesLength > allSubjectFeaturesLength:
            self._normaliseByQuery = False

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        The score is a quotient. In the numerator, we have the number of AA
        locations that are in features that are in hashes that match between
        the subject and the query. The denominator is the number of AA
        locations that are in features which are in all hashes in the matching
        regions of the query and subject.
        Leaving the score like this would mean that a short match can have the
        same score as a long match. To account for this, the quotient from
        above is multiplied by the matched fraction of the shorter sequence.

        @param binIndex: The C{int} index of the bin to examine.
        @return: The C{float} score of that bin.
        """
        # Get the features and their offsets which match in subject and query.
        # These will be used to calculate the numerator of the score.
        matchedQueryFeatures, matchedQueryOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'query')

        matchedSubjectFeatures, matchedSubjectOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'subject')

        # Get the extreme offsets in the matched region of query and subject.
        minQueryOffset = min(matchedQueryOffsets, default=0)
        maxQueryOffset = max(matchedQueryOffsets, default=self._queryLen)
        minSubjectOffset = min(matchedSubjectOffsets, default=0)
        maxSubjectOffset = max(matchedSubjectOffsets, default=self._subjectLen)

        # Get all features and their offsets which are present in the subject
        # and the query within the matched region. These will be used to
        # calculate the denominator.
        unmatchedQueryOffsets = set()
        for feature in filter(
                lambda f: featureInRange(f, minQueryOffset, maxQueryOffset),
                self._allQueryFeatures - matchedQueryFeatures):
            unmatchedQueryOffsets.update(feature.coveredOffsets())
        # The unmatched offsets shouldn't contain any offsets that were
        # matched. This can occur if an unmatched feature overlaps with a
        # matched feature.
        unmatchedQueryOffsets -= matchedQueryOffsets

        unmatchedSubjectOffsets = set()
        for feature in filter(
                lambda f: featureInRange(f, minSubjectOffset,
                                         maxSubjectOffset),
                self._allSubjectFeatures - matchedSubjectFeatures):
            unmatchedSubjectOffsets.update(feature.coveredOffsets())
        # The unmatched offsets shouldn't contain any offsets that were
        # matched. This can occur if an unmatched feature overlaps with a
        # matched feature.
        unmatchedSubjectOffsets -= matchedSubjectOffsets

        matchedOffsetCount = (
            len(matchedQueryOffsets) + len(matchedSubjectOffsets))

        totalOffsetCount = matchedOffsetCount + (
            len(unmatchedQueryOffsets) + len(unmatchedSubjectOffsets))

        # The calculation of the fraction to normalise by length consists of
        # three parts: the numerator is the matchedOffsetCount + either the
        # unmatchedQueryOffsets or the unmatchedSubjectOffsets. The denominator
        # is the numerator + the length of hashes in either subject or query
        # which are outside the matched region. The sequence with less covered
        # indices is used to do the normalisation.

        offsetsNotInMatch = set()
        if self._normaliseByQuery:
            for feature in filterfalse(
                    lambda f: featureInRange(f, minQueryOffset,
                                             maxQueryOffset),
                    self._allQueryFeatures - matchedQueryFeatures):
                offsetsNotInMatch.update(feature.coveredOffsets())
            numerator = matchedOffsetCount + len(unmatchedQueryOffsets)
        else:
            for feature in filterfalse(
                    lambda f: featureInRange(f, minSubjectOffset,
                                             maxSubjectOffset),
                    self._allSubjectFeatures - matchedSubjectFeatures):
                offsetsNotInMatch.update(feature.coveredOffsets())
            numerator = matchedOffsetCount + len(unmatchedSubjectOffsets)

        denominator = numerator + len(offsetsNotInMatch)

        try:
            normaliser = numerator / denominator
        except ZeroDivisionError:
            normaliser = 1.0

        try:
            score = matchedOffsetCount / totalOffsetCount
        except ZeroDivisionError:
            score = 0.0

        return score * normaliser


ALL_SCORE_CLASSES = (MinHashesScore, FeatureMatchingScore, FeatureAAScore)
