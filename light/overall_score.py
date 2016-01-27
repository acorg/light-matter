from light.utils import maxWithDefault, minWithDefault
from light.string import MultilineString
from light.bin_score import (getHashFeatures, featureInRange,
                             histogramBinFeatures)


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

        @return: a C{float} score for the best significant bin (or C{None} if
            there are no significant bins) and a C{dict} with information about
            the score.
        """
        return self._score, self._analysis

    @staticmethod
    def printAnalysis(analysis, margin=''):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
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


def featureIndicesOutsideOffsets(feature, offsets):
    """
    Returns a C{set} containing the offsets of a feature which aren't in the
    offsets specified.

    @param feature: A C{light.features._Feature} subclass (i.e., a
        landmark or a trig point).
    @param offsets: A C{set} of offsets that are covered by a bin.
    @return: A C{set} of featureOffsets that are outside the offsets.
    """
    return feature.coveredOffsets() - offsets


class SignificantBinScore(object):
    """
    Calculate an overall score based on all significant bins. The score is
    calculated as follows:
    The score is a quotient. The numerator consists of the sum of all unique
    offsets in features in pairs that match between subject and query, in all
    significant bins. The denominator consists of the sum of all unique offsets
    in features in pairs in subject and query that don't match, for all
    significant bins.
    The score is corrected by length, by multiplying with another quotient.
    The quotient is calculated for both the subject and the query, and the
    bigger one is chosen. The numerator consists of the offsets of all
    features in the matched region. The denominator consists of all offsets in
    pairs for the whole sequence.

    @param histogram: A C{light.histogram} instance.
    @param significantBins: A C{list} of C{dict}'s where each dict contains
        information about the score, bin and index of a significant bin. This
        list is already sorted by score.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param params: A C{Parameters} instance.
    """
    def __init__(self, histogram, significantBins, query, subject, params):
        self._histogram = histogram
        self._queryLen = len(query)
        self._subjectLen = len(subject)
        self._significantBinIndices = [signBin['index'] for signBin in
                                       significantBins]

        from light.backend import Backend
        backend = Backend()
        backend.configure(params)

        scannedQuery = backend.scan(query)
        allQueryHashes = backend.getHashes(scannedQuery)
        self._allQueryFeatures = getHashFeatures(allQueryHashes)

        scannedSubject = backend.scan(subject)
        allSubjectHashes = backend.getHashes(scannedSubject)
        self._allSubjectFeatures = getHashFeatures(allSubjectHashes)

    def calculateScore(self):
        """
        Calculates the overall score for all significant bins, as described
        above.

        @return: a C{float} score for the best significant bin (or C{None} if
            there are no significant bins) and a C{dict} with information about
            the score.
        """
        if not self._significantBinIndices:
            return None, {}

        # Calculate the first quotient.
        overallMatchedQueryOffsets = set()
        overallMatchedSubjectOffsets = set()

        overallUnMatchedQueryOffsets = set()
        overallUnMatchedSubjectOffsets = set()

        overallQueryStartEnd = []
        overallSubjectStartEnd = []

        # Get the features and their offsets which are matched and unmatched in
        # subject and query in all bins. These will be used to calculate the
        # numerator of the score.
        for binIndex in self._significantBinIndices:
            bin_ = self._histogram[binIndex]
            matchedQueryFeatures, matchedQueryOffsets = histogramBinFeatures(
                bin_, 'query')

            matchedSbjctFeatures, matchedSubjectOffsets = histogramBinFeatures(
                bin_, 'subject')

            # Get the extreme offsets in the matched region of query and
            # subject.
            minQueryOffset = minWithDefault(matchedQueryOffsets,
                                            default=None)
            maxQueryOffset = maxWithDefault(matchedQueryOffsets,
                                            default=None)
            minSubjectOffset = minWithDefault(matchedSubjectOffsets,
                                              default=None)
            maxSubjectOffset = maxWithDefault(matchedSubjectOffsets,
                                              default=None)
            overallQueryStartEnd.append([minQueryOffset, maxQueryOffset])
            overallSubjectStartEnd.append([minSubjectOffset, maxSubjectOffset])

            # Get all features and their offsets which are present in the
            # subject and the query within the matched region. These will be
            # used to calculate the denominator.
            unmatchedQueryOffsets = set()
            for feature in filter(
                    lambda f: featureInRange(f, minQueryOffset,
                                             maxQueryOffset),
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
                    self._allSubjectFeatures - matchedSbjctFeatures):
                unmatchedSubjectOffsets.update(feature.coveredOffsets())
            # The unmatched offsets shouldn't contain any offsets that were
            # matched. This can occur if an unmatched feature overlaps with a
            # matched feature.
            unmatchedSubjectOffsets -= matchedSubjectOffsets

            overallMatchedQueryOffsets.update(matchedQueryOffsets)
            overallMatchedSubjectOffsets.update(matchedSubjectOffsets)
            overallUnMatchedQueryOffsets.update(unmatchedQueryOffsets)
            overallUnMatchedSubjectOffsets.update(unmatchedSubjectOffsets)

        # Calculate the first quotient over all bins.
        matchedOffsetCount = (len(overallMatchedQueryOffsets) +
                              len(overallMatchedSubjectOffsets))
        totalOffsetCount = (matchedOffsetCount +
                            len(overallUnMatchedQueryOffsets) +
                            len(overallUnMatchedSubjectOffsets))
        try:
            matchedRegionScore = matchedOffsetCount / totalOffsetCount
        except ZeroDivisionError:
            matchedRegionScore = 0.0

        # Calculate the score for the second quotient (the length normalizer)
        # Get all offsets that are part of a bin.
        coveredQueryOffsets = set()
        coveredSubjectOffsets = set()
        for start, end in overallQueryStartEnd:
            coveredQueryOffsets.update(range(start, end + 1))
        for start, end in overallSubjectStartEnd:
            coveredSubjectOffsets.update(range(start, end + 1))

        # Get all feature offsets that are not in a bin
        overallOffsetsNotInMatchQuery = set()
        for feature in self._allQueryFeatures:
            coveredFeatIndices = featureIndicesOutsideOffsets(
                feature, coveredQueryOffsets)
            overallOffsetsNotInMatchQuery.update(coveredFeatIndices)
        overallOffsetsNotInMatchQuery -= set(overallMatchedQueryOffsets)
        numeratorQuery = (len(overallMatchedQueryOffsets) +
                          len(overallUnMatchedQueryOffsets))
        denominatorQuery = numeratorQuery + len(overallOffsetsNotInMatchQuery)

        overallOffsetsNotInMatchSubject = set()
        for feature in self._allSubjectFeatures:
            coveredFeatIndices = featureIndicesOutsideOffsets(
                feature, coveredQueryOffsets)
            overallOffsetsNotInMatchSubject.update(coveredFeatIndices)
        overallOffsetsNotInMatchSubject -= set(overallMatchedSubjectOffsets)
        numeratorSubject = (len(overallMatchedSubjectOffsets) +
                            len(overallUnMatchedSubjectOffsets))
        denominatorSubject = (numeratorSubject +
                              len(overallOffsetsNotInMatchSubject))

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
            'matchedSubjectOffsetCount': len(overallMatchedSubjectOffsets),
            'matchedQueryOffsetCount': len(overallMatchedQueryOffsets),
            'matchedRegionScore': matchedRegionScore,
            'numeratorQuery': numeratorQuery,
            'numeratorSubject': numeratorSubject,
            'normaliserQuery': normaliserQuery,
            'normaliserSubject': normaliserSubject,
            'score': score,
            'scoreClass': self.__class__,
            'totalOffsetCount': totalOffsetCount,
            'queryOffsetsInBins': len(coveredQueryOffsets),
            'subjectOffsetsInBins': len(coveredSubjectOffsets),
        }

        return score, analysis

    @staticmethod
    def printAnalysis(analysis, margin=''):
        """
        Convert an analysis to a nicely formatted string.

        @param analysis: A C{dict} with information about the score and its
            calculation.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} human-readable version of the last analysis.
        """
        result = MultilineString(margin=margin)

        result.extend([
            'Overall score method: %s' % analysis['scoreClass'].__name__,
            'Overall score: %s' % analysis['score'],
            ('Total (query+subject) AA offsets in matched pairs in all bins: '
             '%(matchedOffsetCount)d' % analysis),
            ('Subject AA offsets in matched pairs in all bins: '
             '%(matchedSubjectOffsetCount)d' % analysis),
            ('Query AA offsets in matched pairs in all bins: '
             '%(matchedQueryOffsetCount)d' % analysis),
            ('Total (query+subject) AA offsets in hashes in matched region: '
             '%(totalOffsetCount)d' % analysis),
            ('Matched region score %(matchedRegionScore).4f '
             '(%(matchedOffsetCount)d / %(totalOffsetCount)d)' % analysis),
            ('Query normalizer: %(normaliserQuery).4f (%(numeratorQuery)d / '
             '%(denominatorQuery)d)' % analysis),
            ('Subject normalizer: %(normaliserSubject).4f '
             '(%(numeratorSubject)d / %(denominatorSubject)d)' % analysis),
            ('Total query offsets that are in a bin: '
             '%(queryOffsetsInBins)d' % analysis),
            ('Total subject offsets that are in a bin: '
             '%(subjectOffsetsInBins)d' % analysis),
        ])

        return str(result)

ALL_OVERALL_SCORE_CLASSES = [BestBinScore, SignificantBinScore]
