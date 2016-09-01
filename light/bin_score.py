from __future__ import division

from collections import defaultdict
from six.moves import filterfalse
from warnings import warn

from light.string import MultilineString
from light.utils import maxWithDefault, minWithDefault


class NoneScore(object):
    """
    Unconditionally set bin scores to C{None}.

    This is useful when the individual scores of bins aren't needed, most
    likely because an overall match score between query and subject is
    more useful.
    """
    def __init__(self):
        self._analysis = {
            'score': None,
            'scoreClass': self.__class__,
        }

    def calculateScore(self, binIndex):
        """
        Calculate the score (always C{None}) for a histogram bin.

        @param binIndex: The C{int} index of the bin to (supposedly) examine.
        @return: A 2-tuple, containing C{None} (the score of the bin) and a
            C{dict} with the analysis leading to the score.
        """
        return None, self._analysis


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
        self._minHashCount = minHashCount
        self.analysis = None

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A 2-tuple, containing the C{float} score of the bin and a
            C{dict} with the analysis leading to the score.
        """
        binCount = len(self._histogram[binIndex])
        try:
            score = float(binCount) / self._minHashCount
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
                 (binCount, self._minHashCount), RuntimeWarning)

        analysis = {
            'minHashCount': self._minHashCount,
            'binCount': binCount,
            'score': score,
            'scoreClass': self.__class__,
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin='', result=None):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} human-readable
            version of the last analysis, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.extend([
            'Score method: %s' % analysis['scoreClass'].__name__,
            'Minimum hash count: %(minHashCount)d' % analysis,
            'Bin count: %(binCount)d' % analysis,
            'Score: %(score).4f' % analysis,
        ])

        if not returnNone:
            return str(result)


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


def weightedHistogramBinFeatures(bin_, queryOrSubject, weights):
    """
    Extract all features of type C{queryOrSubject} (a C{str}, either
    'query' or 'subject') from a bin and return them as a set. Also returned is
    a dictionary, where each key is an offset and the values is a list of
    weights associated with that position (it has to be a list, because if
    landmarks overlap, an offset may have two weights).

    @param bin_: A C{light.histogram.Histogram} bin.
    @param queryOrSubject: A C{str}, to indicate which features to extract,
        either 'query' or 'subject'.
    @param weights: A C{dict} specifying the weight that should be give to
        each landmark.
    @raise KeyError: If C{queryOrSubject} is not 'query' or 'subject'.
    @return: A 2-tuple, containing 1) a C{set} of all features (landmarks and
        trig points) in the hashes in C{bin_}, 2) a C{set} of all offsets of
        all features (of type C{queryOrSubject}) in the bin.
    """
    allFeatures = set()
    allOffsets = defaultdict(list)
    # There is no error checking that queryOrSubject is 'query' or
    # 'subject' as the following will raise a KeyError if it cannot access
    # the dict key in the bin element.
    for hashInfo in bin_:
        for suffix in 'Landmark', 'TrigPoint':
            feature = hashInfo[queryOrSubject + suffix]
            weight = weights[feature.name]
            allFeatures.add(feature)
            for offset in feature.coveredOffsets():
                allOffsets[offset].append(weight)

    return allFeatures, allOffsets


def featureInRange(feature, minOffset, maxOffset):
    """
    Does a feature fall within a min/max offset range?

    @param feature: A C{light.features._Feature} subclass (i.e., a
        landmark or a trig point).
    @param minOffset: The minimum allowed offset for the feature start. If
        C{None} is passed, no feature is considered in range.
    @param maxOffset: The maximum allowed offset for the feature end.
    @return: A C{bool} to indicate whether the feature falls completely
        within the allowed range.
    """
    return minOffset is not None and feature.offset >= minOffset and (
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
        for lmTpPair in hashInfo:
            allFeatures.update(lmTpPair)

    return allFeatures


class FeatureMatchingScore(object):
    """
    Calculates the score for histogram bins based on the features present
    within the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param dbParams: A C{DatabaseParameters} instance.
    @param findParams: An instance of C{light.parameters.FindParameters} or
        C{None} to use default find parameters.
    """
    def __init__(self, histogram, query, subject, dbParams, findParams=None):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)
        from light.parameters import FindParameters
        self._findParams = findParams or FindParameters()
        from light.backend import Backend
        backend = Backend()
        backend.configure(dbParams)
        scannedQuery = backend.scan(query)
        self._allQueryFeatures = set(scannedQuery.landmarks +
                                     scannedQuery.trigPoints)
        scannedSubject = backend.scan(subject.read)
        self._allSubjectFeatures = set(scannedSubject.landmarks +
                                       scannedSubject.trigPoints)

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A 2-tuple, containing the C{float} score of the bin and a
            C{dict} with the analysis leading to the score.
        """
        queryFeatures, queryOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'query')
        subjectFeatures, subjectOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'subject')

        matchScore = self._findParams.featureMatchScore * (
            len(queryFeatures) + len(subjectFeatures))

        minQueryOffset = minWithDefault(queryOffsets, default=None)
        maxQueryOffset = maxWithDefault(queryOffsets, default=None)
        minSubjectOffset = minWithDefault(subjectOffsets, default=None)
        maxSubjectOffset = maxWithDefault(subjectOffsets, default=None)

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

        score = matchScore + mismatchScore

        # We could put more in here, but I have a feeling we wont be using
        # this score method as FeatureAAScore is (hopefully) better. If we
        # need to add more detail (e.g., the number of features in each
        # class) we can easily add it.
        analysis = {
            'minQueryOffset': minQueryOffset,
            'maxQueryOffset': maxQueryOffset,
            'minSubjectOffset': minSubjectOffset,
            'maxSubjectOffset': maxSubjectOffset,
            'matchScore': matchScore,
            'mismatchScore': mismatchScore,
            'score': score,
            'scoreClass': self.__class__,
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin='', result=None):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} human-readable
            version of the last analysis, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.extend([
            'Score method: %s' % analysis['scoreClass'].__name__,
            ('Matched offset range in query: %(minQueryOffset)d to '
             '%(maxQueryOffset)d' % analysis),
            ('Matched offset range in subject: %(minSubjectOffset)d to '
             '%(maxSubjectOffset)d' % analysis),
            ('Match score: %(matchScore).4f' % analysis),
            ('Mismatch score: %(mismatchScore).4f' % analysis),
            'Score: %(score).4f' % analysis,
        ])

        if not returnNone:
            return str(result)


class FeatureAAScore(object):
    """
    Calculates the score for histogram bins based on the count of amino acids
    present in the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param dbParams: A C{DatabaseParameters} instance.
    """
    def __init__(self, histogram, query, subject, dbParams):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)

        from light.backend import Backend
        backend = Backend()
        backend.configure(dbParams)

        scannedQuery = backend.scan(query)
        allQueryHashes = backend.getHashes(scannedQuery)
        self._allQueryFeatures = getHashFeatures(allQueryHashes)

        scannedSubject = backend.scan(subject.read)
        allSubjectHashes = backend.getHashes(scannedSubject)
        self._allSubjectFeatures = getHashFeatures(allSubjectHashes)

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        The score is the product of two quotients. The first quotient is the
        matched region score (MRS). In the numerator, we have the number of AA
        locations that are in features that are in hashes that match between
        the subject and the query. The denominator is the number of AA
        locations that are in features which are in all hashes in the matching
        regions of the query and subject.
        Leaving the score like this would mean that a short match can have the
        same score as a long match. To account for this we multiply the MRS
        with a length normaliser (LN). The LN is the quotient of all AA
        locations in hashes in the matched region in the subject or the query
        divided by all AA locations in hashes in the subject or the query,
        whichever of the two is bigger.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A 2-tuple, containing the C{float} score of the bin and a
            C{dict} with the analysis leading to the score.
        """
        # Get the features and their offsets which match in subject and query.
        # These will be used to calculate the numerator of the score.
        matchedQueryFeatures, matchedQueryOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'query')

        matchedSubjectFeatures, matchedSubjectOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'subject')

        # Get the extreme offsets in the matched region of query and subject.
        minQueryOffset = minWithDefault(matchedQueryOffsets, default=None)
        maxQueryOffset = maxWithDefault(matchedQueryOffsets, default=None)
        minSubjectOffset = minWithDefault(matchedSubjectOffsets, default=None)
        maxSubjectOffset = maxWithDefault(matchedSubjectOffsets, default=None)

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

        try:
            matchedRegionScore = matchedOffsetCount / totalOffsetCount
        except ZeroDivisionError:
            matchedRegionScore = 0.0

        # The calculation of the fraction to normalise by length consists of
        # three parts: the numerator is the matchedOffsetCount + either the
        # unmatchedQueryOffsets or the unmatchedSubjectOffsets. The denominator
        # is the numerator + the length of hashes in either subject or query
        # which are outside the matched region. The sequence with less covered
        # indices is used to do the normalisation.

        offsetsNotInMatchQuery = set()
        for feature in filterfalse(
                lambda f: featureInRange(f, minQueryOffset,
                                         maxQueryOffset),
                self._allQueryFeatures - matchedQueryFeatures):
            offsetsNotInMatchQuery.update(feature.coveredOffsets())
        offsetsNotInMatchQuery -= matchedQueryOffsets
        numeratorQuery = (len(matchedQueryOffsets) +
                          len(unmatchedQueryOffsets))
        denominatorQuery = numeratorQuery + len(offsetsNotInMatchQuery)

        offsetsNotInMatchSubject = set()
        for feature in filterfalse(
                lambda f: featureInRange(f, minSubjectOffset,
                                         maxSubjectOffset),
                self._allSubjectFeatures - matchedSubjectFeatures):
            offsetsNotInMatchSubject.update(feature.coveredOffsets())
        offsetsNotInMatchSubject -= matchedSubjectOffsets
        numeratorSubject = (len(matchedSubjectOffsets) +
                            len(unmatchedSubjectOffsets))
        denominatorSubject = numeratorSubject + len(offsetsNotInMatchSubject)

        try:
            normaliserQuery = numeratorQuery / denominatorQuery
        except ZeroDivisionError:
            normaliserQuery = 1.0

        try:
            normaliserSubject = numeratorSubject / denominatorSubject
        except ZeroDivisionError:
            normaliserSubject = 1.0

        score = matchedRegionScore * max(normaliserQuery, normaliserSubject)

        analysis = {
            'denominatorQuery': denominatorQuery,
            'denominatorSubject': denominatorSubject,
            'matchedOffsetCount': matchedOffsetCount,
            'matchedSubjectOffsetCount': len(matchedSubjectOffsets),
            'matchedQueryOffsetCount': len(matchedQueryOffsets),
            'matchedRegionScore': matchedRegionScore,
            'maxQueryOffset': maxQueryOffset,
            'maxSubjectOffset': maxSubjectOffset,
            'minQueryOffset': minQueryOffset,
            'minSubjectOffset': minSubjectOffset,
            'numeratorQuery': numeratorQuery,
            'numeratorSubject': numeratorSubject,
            'normaliserQuery': normaliserQuery,
            'normaliserSubject': normaliserSubject,
            'score': score,
            'scoreClass': self.__class__,
            'totalOffsetCount': totalOffsetCount,
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin='', result=None):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} human-readable
            version of the last analysis, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.extend([
            'Score method: %s' % analysis['scoreClass'].__name__,
            ('Matched offset range in query: %(minQueryOffset)d to '
             '%(maxQueryOffset)d' % analysis),
            ('Matched offset range in subject: %(minSubjectOffset)d to '
             '%(maxSubjectOffset)d' % analysis),
            ('Total (query+subject) AA offsets in matched hashes: '
             '%(matchedOffsetCount)d' % analysis),
            ('Subject AA offsets in matched hashes: '
             '%(matchedSubjectOffsetCount)d' % analysis),
            ('Query AA offsets in matched hashes: '
             '%(matchedQueryOffsetCount)d' % analysis),
            ('Total (query+subject) AA offsets in hashes in matched region: '
             '%(totalOffsetCount)d' % analysis),
            ('Matched region score %(matchedRegionScore).4f '
             '(%(matchedOffsetCount)d / %(totalOffsetCount)d)' % analysis),
            ('Query normalizer: %(normaliserQuery).4f (%(numeratorQuery)d / '
             '%(denominatorQuery)d)' % analysis),
            ('Subject normalizer: %(normaliserSubject).4f '
             '(%(numeratorSubject)d / %(denominatorSubject)d)' % analysis),
            'Score: %(score).4f' % analysis,
        ])

        if returnNone is not None:
            return str(result)


def getWeightedOffsets(offsetDict):
    """
    Calculate the weighted offsets.

    @param offsetDict: A C{dict} where the key is an offset and the value is
        the weight associated with that offset.
    @return: A C{float} weighted offset count.
    """
    return sum(max(weights) for weights in offsetDict.values())


class WeightedFeatureAAScore(object):
    """
    Calculates the score for histogram bins based on the count of amino acids
    present in the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).
    Weight the score by importance of the features in the matched region.

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param dbParams: A C{DatabaseParameters} instance.
    @param weights: If not C{None}, a C{dict} of weights to be given to each
        feature.
    """
    def __init__(self, histogram, query, subject, dbParams, weights=None):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)

        self._weights = self.DEFAULT_WEIGHTS if weights is None else weights

        from light.backend import Backend
        backend = Backend()
        backend.configure(dbParams)

        scannedQuery = backend.scan(query)
        allQueryHashes = backend.getHashes(scannedQuery)
        self._allQueryFeatures = getHashFeatures(allQueryHashes)

        scannedSubject = backend.scan(subject.read)
        allSubjectHashes = backend.getHashes(scannedSubject)
        self._allSubjectFeatures = getHashFeatures(allSubjectHashes)

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        The score is a quotient. In the numerator, we have the weighted number
        of AA locations that are in features that are in hashes that match
        between the subject and the query. The denominator is the weighted
        number of AA locations that are in features which are in all hashes in
        the matching regions of the query and subject.
        Leaving the score like this would mean that a short match can have the
        same score as a long match. To account for this, the quotient from
        above is multiplied by the matched fraction of the shorter sequence.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A 2-tuple, containing the C{float} score of the bin and a
            C{dict} with the analysis leading to the score.
        """
        # Get the features and their offsets with associated weights which
        # match in subject and query.
        # These will be used to calculate the numerator of the score.
        matchedQFeatures, matchedQOffsets = weightedHistogramBinFeatures(
            self._histogram[binIndex], 'query', self._weights)

        matchedSFeatures, matchedSOffsets = weightedHistogramBinFeatures(
            self._histogram[binIndex], 'subject', self._weights)

        # Get the extreme offsets in the matched region of query and subject.
        minQueryOffset = minWithDefault(matchedQOffsets, default=None)
        maxQueryOffset = maxWithDefault(matchedQOffsets, default=None)
        minSubjectOffset = minWithDefault(matchedSOffsets, default=None)
        maxSubjectOffset = maxWithDefault(matchedSOffsets, default=None)

        # Get all features and their offsets with their associated weights
        # which are present in the subject and the query within the matched
        # region.
        # These will be used to calculate the denominator of the score.
        unmatchedQueryOffsets = defaultdict(list)
        for feature in filter(
                lambda f: featureInRange(f, minQueryOffset, maxQueryOffset),
                self._allQueryFeatures - matchedQFeatures):
            for offset in feature.coveredOffsets():
                unmatchedQueryOffsets[offset].append(
                    self._weights[feature.name])

        unmatchedSubjectOffsets = defaultdict(list)
        for feature in filter(
                lambda f: featureInRange(f, minSubjectOffset,
                                         maxSubjectOffset),
                self._allSubjectFeatures - matchedSFeatures):
            for offset in feature.coveredOffsets():
                unmatchedSubjectOffsets[offset].append(
                    self._weights[feature.name])

        # The unmatched offsets in the query and the subject shouldn't contain
        # any offsets that were matched. This can occur if an unmatched feature
        # overlaps with a matched feature.
        for offset in matchedQOffsets.keys():
            unmatchedQueryOffsets.pop(offset, None)

        for offset in matchedSOffsets.keys():
            unmatchedSubjectOffsets.pop(offset, None)

        matchedWeightsCount = (
            getWeightedOffsets(matchedQOffsets) +
            getWeightedOffsets(matchedSOffsets))

        totalWeightsCount = matchedWeightsCount + (
            getWeightedOffsets(unmatchedQueryOffsets) +
            getWeightedOffsets(unmatchedSubjectOffsets))

        # Calculate the weighted score of the features within the matched
        # region.
        try:
            matchedRegionScore = matchedWeightsCount / totalWeightsCount
        except ZeroDivisionError:
            matchedRegionScore = 0.0

        # The calculation of the fraction to normalise by length consists of
        # three parts: the numerator is the matchedOffsetCount + either the
        # unmatchedQueryOffsets or the unmatchedSubjectOffsets. The denominator
        # is the numerator + the length of hashes in either subject or query
        # which are outside the matched region. The sequence with less covered
        # indices is used to do the normalisation.
        # Note that the fraction to normalise by length is not weighted.

        offsetsNotInMatchQuery = set()
        for feature in filterfalse(
                lambda f: featureInRange(f, minQueryOffset,
                                         maxQueryOffset),
                self._allQueryFeatures - matchedQFeatures):
            offsetsNotInMatchQuery.update(feature.coveredOffsets())

        matchedQOffsetsSet = set(matchedQOffsets.keys())
        offsetsNotInMatchQuery -= matchedQOffsetsSet

        numeratorQuery = (len(matchedQOffsets) +
                          len(unmatchedQueryOffsets))
        denominatorQuery = numeratorQuery + len(offsetsNotInMatchQuery)

        offsetsNotInMatchSubject = set()
        for feature in filterfalse(
                lambda f: featureInRange(f, minSubjectOffset,
                                         maxSubjectOffset),
                self._allSubjectFeatures - matchedSFeatures):
            offsetsNotInMatchSubject.update(feature.coveredOffsets())

        matchedSOffsetsSet = set(matchedSOffsets.keys())
        offsetsNotInMatchSubject -= matchedSOffsetsSet

        numeratorSubject = (len(matchedSOffsets) +
                            len(unmatchedSubjectOffsets))
        denominatorSubject = numeratorSubject + len(offsetsNotInMatchSubject)

        # Calculate the fraction to normalise by length.
        try:
            normaliserQuery = numeratorQuery / denominatorQuery
        except ZeroDivisionError:
            normaliserQuery = 1.0

        try:
            normaliserSubject = numeratorSubject / denominatorSubject
        except ZeroDivisionError:
            normaliserSubject = 1.0

        # Calculate the final score.
        score = matchedRegionScore * max(normaliserQuery, normaliserSubject)

        analysis = {
            'denominatorQuery': denominatorQuery,
            'denominatorSubject': denominatorSubject,
            'matchedOffsetCount': matchedWeightsCount,
            'matchedSubjectOffsetCount': len(matchedSOffsets),
            'matchedQueryOffsetCount': len(matchedQOffsets),
            'weightedMatchedQueryOffsetCount':
                getWeightedOffsets(matchedQOffsets),
            'weightedMatchedSubjectOffsetCount':
                getWeightedOffsets(matchedSOffsets),
            'matchedRegionScore': matchedRegionScore,
            'maxQueryOffset': maxQueryOffset,
            'maxSubjectOffset': maxSubjectOffset,
            'minQueryOffset': minQueryOffset,
            'minSubjectOffset': minSubjectOffset,
            'numeratorQuery': numeratorQuery,
            'numeratorSubject': numeratorSubject,
            'normaliserQuery': normaliserQuery,
            'normaliserSubject': normaliserSubject,
            'score': score,
            'scoreClass': self.__class__,
            'totalOffsetCount': totalWeightsCount,
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin='', result=None):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} human-readable
            version of the last analysis, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.extend([
            'Score method: %s' % analysis['scoreClass'].__name__,
            ('Matched offset range in query: %(minQueryOffset)d to '
             '%(maxQueryOffset)d' % analysis),
            ('Matched offset range in subject: %(minSubjectOffset)d to '
             '%(maxSubjectOffset)d' % analysis),
            ('Total (query+subject) AA offsets in matched hashes: '
             '%(matchedOffsetCount)d' % analysis),
            ('Subject AA offsets in matched hashes: '
             '%(matchedSubjectOffsetCount)d' % analysis),
            ('Query AA offsets in matched hashes: '
             '%(matchedQueryOffsetCount)d' % analysis),
            ('Total (query+subject) AA offsets in hashes in matched region: '
             '%(totalOffsetCount)d' % analysis),
            ('Weighted Subject AA offsets in matched hashes: '
             '%(weightedMatchedSubjectOffsetCount)d' % analysis),
            ('Weighted Query AA offsets in matched hashes: '
             '%(weightedMatchedQueryOffsetCount)d' % analysis),
            ('Matched region score %(matchedRegionScore).4f '
             '(%(matchedOffsetCount)d / %(totalOffsetCount)d)' % analysis),
            ('Query normalizer: %(normaliserQuery).4f (%(numeratorQuery)d / '
             '%(denominatorQuery)d)' % analysis),
            ('Subject normalizer: %(normaliserSubject).4f '
             '(%(numeratorSubject)d / %(denominatorSubject)d)' % analysis),
            'Score: %(score).4f' % analysis,
        ])

        if not returnNone:
            return str(result)


class FeatureAALengthScore(object):
    """
    Calculates the score for histogram bins based on the count of amino acids
    present in the regions of the query and subject that had a significant
    alignment (i.e., caused a histogram bin to be considered significant).
    Note that the FeatureAALengthScore and the FeatureAAScore are very similar.
    The only difference is in the calculation of the length normaliser (see
    below).

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param dbParams: A C{DatabaseParameters} instance.
    """
    def __init__(self, histogram, query, subject, dbParams):
        self._histogram = histogram
        self.queryLen = len(query)
        self.subjectLen = len(subject)

        from light.backend import Backend
        backend = Backend()
        backend.configure(dbParams)

        scannedQuery = backend.scan(query)
        allQueryHashes = backend.getHashes(scannedQuery)
        self._allQueryFeatures = getHashFeatures(allQueryHashes)

        scannedSubject = backend.scan(subject.read)
        allSubjectHashes = backend.getHashes(scannedSubject)
        self._allSubjectFeatures = getHashFeatures(allSubjectHashes)

    def calculateScore(self, binIndex):
        """
        Calculates the score for a given histogram bin.

        The score is the product of two quotients. The first quotient is the
        matched region score (MRS). In the numerator, we have the number of AA
        locations that are in features that are in hashes that match between
        the subject and the query. The denominator is the number of AA
        locations that are in features which are in all hashes in the matching
        regions of the query and subject.
        Leaving the score like this would mean that a short match can have the
        same score as a long match. To account for this we multiply the MRS
        with a length normaliser (LN). The LN is the quotient of all AA
        locations in the matched region divided by the total length of the
        subject or the query, whichever of the two is smaller.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A 2-tuple, containing the C{float} score of the bin and a
            C{dict} with the analysis leading to the score.
        """
        # Get the features and their offsets which match in subject and query.
        # These will be used to calculate the numerator of the score.
        matchedQueryFeatures, matchedQueryOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'query')

        matchedSubjectFeatures, matchedSubjectOffsets = histogramBinFeatures(
            self._histogram[binIndex], 'subject')

        # Get the extreme offsets in the matched region of query and subject.
        minQueryOffset = minWithDefault(matchedQueryOffsets, default=None)
        maxQueryOffset = maxWithDefault(matchedQueryOffsets, default=None)
        minSubjectOffset = minWithDefault(matchedSubjectOffsets, default=None)
        maxSubjectOffset = maxWithDefault(matchedSubjectOffsets, default=None)

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
        try:
            matchedRegionScore = matchedOffsetCount / totalOffsetCount
        except ZeroDivisionError:
            matchedRegionScore = 0.0

        # The length normaliser is a quotient. The numerator is the number of
        # AA in the matched region. The denominator is the length of the
        # subject or the query, whichever one is smaller.

        if maxQueryOffset:
            numerator = maxQueryOffset - minQueryOffset + 1
        else:
            numerator = 0

        lengthNormaliser = numerator / min(self.queryLen, self.subjectLen)

        score = matchedRegionScore * lengthNormaliser

        analysis = {
            'denominatorQuery': self.queryLen,
            'denominatorSubject': self.subjectLen,
            'matchedOffsetCount': matchedOffsetCount,
            'matchedSubjectOffsetCount': len(matchedSubjectOffsets),
            'matchedQueryOffsetCount': len(matchedQueryOffsets),
            'matchedRegionScore': matchedRegionScore,
            'maxQueryOffset': maxQueryOffset,
            'maxSubjectOffset': maxSubjectOffset,
            'minQueryOffset': minQueryOffset,
            'minSubjectOffset': minSubjectOffset,
            'matchedRegionSize': numerator,
            'normaliserQuery': numerator / self.queryLen,
            'normaliserSubject': numerator / self.subjectLen,
            'score': score,
            'scoreClass': self.__class__,
            'totalOffsetCount': totalOffsetCount,
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin='', result=None):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} human-readable
            version of the last analysis, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.extend([
            'Score method: %s' % analysis['scoreClass'].__name__,
            ('Matched offset range in query: %(minQueryOffset)d to '
             '%(maxQueryOffset)d' % analysis),
            ('Matched offset range in subject: %(minSubjectOffset)d to '
             '%(maxSubjectOffset)d' % analysis),
            ('Total (query+subject) AA offsets in matched hashes: '
             '%(matchedOffsetCount)d' % analysis),
            ('Subject AA offsets in matched hashes: '
             '%(matchedSubjectOffsetCount)d' % analysis),
            ('Query AA offsets in matched hashes: '
             '%(matchedQueryOffsetCount)d' % analysis),
            ('Total (query+subject) AA offsets in hashes in matched region: '
             '%(totalOffsetCount)d' % analysis),
            ('Matched region score %(matchedRegionScore).4f '
             '(%(matchedOffsetCount)d / %(totalOffsetCount)d)' % analysis),
            ('Query normalizer: %(normaliserQuery).4f (%(matchedRegionSize)d '
             '/ %(denominatorQuery)d)' % analysis),
            ('Subject normalizer: %(normaliserSubject).4f '
             '(%(matchedRegionSize)d / %(denominatorSubject)d)' % analysis),
            'Score: %(score).4f' % analysis,
        ])

        if returnNone is not None:
            return str(result)

ALL_BIN_SCORE_CLASSES = (NoneScore, MinHashesScore, FeatureMatchingScore,
                         FeatureAAScore, WeightedFeatureAAScore,
                         FeatureAALengthScore)
