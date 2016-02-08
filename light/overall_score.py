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


def offsetsInBin(bin_, queryOrSubject, allFeatures):
    """
    Calculate sets of matched and unmatched offsets inside the bin, as well as
    the minimum and maximum of of the matched offsets.

    @param bin_: A C{light.histogram.Histogram} bin.
    @param queryOrSubject: A C{str}, to indicate which features to extract,
        either 'query' or 'subject'.
     @param allFeatures: All the features that were found in the sequence
        in question (either the query or the subject).
    @raise KeyError: If C{queryOrSubject} is not 'query' or 'subject'.
    @return: A 4-tuple, with 1) the set of offsets of features that were
        matched in the bin, 2) the set of offsets in features in match
        region which were not in the match, 3) The minimum of the
        matched offsets, 4) the maximum of the matched offsets.
    """
    features, offsets = histogramBinFeatures(bin_, queryOrSubject)

    # Get the min/max offsets of features in the matched region.
    minOffset = min(offsets)
    maxOffset = max(offsets)

    # Get the offsets of all features which are within the matched
    # region but which are not matched in this bin.
    unmatchedOffsets = set()
    for feature in filter(
            lambda f: featureInRange(f, minOffset, maxOffset),
            allFeatures - features):
        unmatchedOffsets.update(feature.coveredOffsets())

    # The unmatched offsets shouldn't contain any offsets that were
    # matched. (This occur when an unmatched feature overlaps with
    # a matched feature.)
    unmatchedOffsets -= offsets

    return offsets, unmatchedOffsets, minOffset, maxOffset


def computeLengthNormalizer(allFeatures, overallMatchedOffsets,
                            overallUnmatchedOffsets, offsetsInBins):
    """
    Compute a Length Normalizer (LN) for the SignificantBinScore class.

    @param allFeatures: All the features that were found in the sequence
        in question (either the query or the subject).
    @param overallMatchedOffsets: A C{set} of all the offsets in features
        that were in the match.
    @param overallUnmatchedoffsets: A C{set} of all the offsets in features
        that were in the matched region but which were not part of the match.
    @param offsetsInBins: A C{set} of all the offsets in features across all
        significant bins.
    @return: A 3-tuple, containing:
        1) the C{float} length normalizer in the range [0.0, 1.0],
        2) the C{int} numerator of the calculation, and
        3) the C{int} denominator of the calculation.
    """
    # Get all feature offsets that are not in any bin.
    overallOffsetsNotInMatches = set()
    update = overallOffsetsNotInMatches.update
    for feature in allFeatures:
        update(feature.coveredOffsets() - offsetsInBins)

    numerator = len(overallMatchedOffsets) + len(overallUnmatchedOffsets)
    denominator = numerator + len(overallOffsetsNotInMatches)

    try:
        normalizer = numerator / denominator
    except ZeroDivisionError:
        normalizer = 1.0

    return normalizer, numerator, denominator


class SignificantBinScore(object):
    """
    Calculate an overall score based on all significant bins.

    The overall score is calculated as a product:

        score = MRS * LN

        where MRS is a C{float} Matched Region Score, and
        LN is a C{float} Length Normalizer.

    The overall score is always a C{float} in the range [0.0 to 1.0].

    MRS is a quotient. The numerator is the number of all unique offsets in
    features in pairs that match between subject and query, in all
    significant bins. The denominator is the number of all unique offsets
    in features in pairs in subject and query that don't match, for all
    significant bins plus the numerator.

    MRS is normalized by length, by multiplying it by LN, another quotient
    derived from either the query or subject.  A quotient is calculated
    individually for the subject and the query, and the larger is used as
    LN. The numerator is the number of unique offsets in all features in the
    matched region in the query (subject).  The denominator is the number of
    unique offsets in all features in the whole of the query (subject)
    sequence.

    @param significantBins: A C{list} of C{dict}s where each C{dict} contains
        information about the score, bin and index of a significant bin. This
        list is already sorted by bin score.
    @param query: A C{dark.reads.AARead} instance.
    @param subject: A C{light.subject.Subject} instance (a subclass of
        C{dark.reads.AARead}).
    @param dbParams: A C{DatabaseParameters} instance.
    """
    def __init__(self, significantBins, query, subject, dbParams):
        # The only thing we do here is store what we have been passed. All
        # work happens in calculateScore (which is only called once).
        self._significantBins = significantBins
        self._query = query
        self._subject = subject
        self._dbParams = dbParams

    def calculateScore(self):
        """
        Calculates the overall score for all significant bins, as described
        above.

        @return: a C{float} overall score for all significant bins (or C{None}
            if there are no significant bins) and a C{dict} with information
            about the score.
        """
        if not self._significantBins:
            analysis = {
                'score': None,
                'scoreClass': self.__class__,
            }

            return None, analysis

        from light.backend import Backend
        backend = Backend()
        backend.configure(self._dbParams)

        allQueryFeatures = getHashFeatures(backend.getHashes(
            backend.scan(self._query)))

        allSubjectFeatures = getHashFeatures(backend.getHashes(
            backend.scan(self._subject)))

        # overallMatchedQueryOffsets and overallMatchedSubjectOffsets will
        # contain all int offsets that are in matching features (and thus
        # inside the matched region).
        overallMatchedQueryOffsets = set()
        overallMatchedSubjectOffsets = set()

        # overallUnmatchedQueryOffsets and overallUnmatchedSubjectOffsets
        # will contain all int offsets that are in features that don't match,
        # but which are inside the matched region.
        overallUnmatchedQueryOffsets = set()
        overallUnmatchedSubjectOffsets = set()

        # The set of all offsets in all bins (whether or not the offsets are in
        # matched features, unmatched features, or not in any feature.
        queryOffsetsInBins = set()
        subjectOffsetsInBins = set()

        # Get the features and their offsets which are matched and unmatched in
        # subject and query in all bins.
        for bin_ in (sb['bin'] for sb in self._significantBins):
            # Query.
            matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
                offsetsInBin(bin_, 'query', allQueryFeatures))
            overallMatchedQueryOffsets.update(matchedOffsets)
            overallUnmatchedQueryOffsets.update(unmatchedOffsets)
            queryOffsetsInBins.update(range(minOffset, maxOffset + 1))

            # Subject.
            matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
                offsetsInBin(bin_, 'subject', allSubjectFeatures))
            overallMatchedSubjectOffsets.update(matchedOffsets)
            overallUnmatchedSubjectOffsets.update(unmatchedOffsets)
            subjectOffsetsInBins.update(range(minOffset, maxOffset + 1))

        # Make sure none of the overall matched offsets are in the overall
        # unmatchedOffsets.
        overallMatchedQueryOffsets -= overallUnmatchedQueryOffsets
        overallMatchedSubjectOffsets -= overallUnmatchedSubjectOffsets

        # Overall score calculation step 1: the matched region score (MRS).
        matchedOffsetCount = (len(overallMatchedQueryOffsets) +
                              len(overallMatchedSubjectOffsets))
        totalOffsetCount = (matchedOffsetCount +
                            len(overallUnmatchedQueryOffsets) +
                            len(overallUnmatchedSubjectOffsets))

        try:
            matchedRegionScore = matchedOffsetCount / totalOffsetCount
        except ZeroDivisionError:
            # A small optimization could be done here. If the MRS is zero,
            # we already know the overall score will be zero, so we could
            # return at this point. To keep things simple, for now, just
            # continue with the overall calculation.
            matchedRegionScore = 0.0

        # Overall score calculation step 2: the length normalizer (LN).

        normalizerQuery, numeratorQuery, denominatorQuery = (
            computeLengthNormalizer(
                allQueryFeatures, overallMatchedQueryOffsets,
                overallUnmatchedQueryOffsets, queryOffsetsInBins))

        # There is a small optimization that could be done at this point.
        # If the query normalizer is 1.0, don't bother to compute a
        # normalizer for the subject (due to the use of max() below and
        # because a normalizer is always <= 1.0).  But to keep the code
        # simpler, for now, we still compute both normalizers.

        normalizerSubject, numeratorSubject, denominatorSubject = (
            computeLengthNormalizer(
                allSubjectFeatures, overallMatchedSubjectOffsets,
                overallUnmatchedSubjectOffsets, subjectOffsetsInBins))

        # Calculate the final score, as descibed in the docstring.
        score = matchedRegionScore * max(normalizerQuery, normalizerSubject)

        analysis = {
            'denominatorQuery': denominatorQuery,
            'denominatorSubject': denominatorSubject,
            'matchedOffsetCount': matchedOffsetCount,
            'matchedSubjectOffsetCount': len(overallMatchedSubjectOffsets),
            'matchedQueryOffsetCount': len(overallMatchedQueryOffsets),
            'matchedRegionScore': matchedRegionScore,
            'numeratorQuery': numeratorQuery,
            'numeratorSubject': numeratorSubject,
            'normalizerQuery': normalizerQuery,
            'normalizerSubject': normalizerSubject,
            'score': score,
            'scoreClass': self.__class__,
            'totalOffsetCount': totalOffsetCount,
            'queryOffsetsInBins': len(queryOffsetsInBins),
            'subjectOffsetsInBins': len(subjectOffsetsInBins),
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
            ('Query normalizer: %(normalizerQuery).4f (%(numeratorQuery)d / '
             '%(denominatorQuery)d)' % analysis),
            ('Subject normalizer: %(normalizerSubject).4f '
             '(%(numeratorSubject)d / %(denominatorSubject)d)' % analysis),
            ('Total query offsets that are in a bin: '
             '%(queryOffsetsInBins)d' % analysis),
            ('Total subject offsets that are in a bin: '
             '%(subjectOffsetsInBins)d' % analysis),
        ])

        return str(result)

ALL_OVERALL_SCORE_CLASSES = [BestBinScore, SignificantBinScore]
