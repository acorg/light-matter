from __future__ import division

import six
from dark.aa import PROPERTY_CLUSTERS

from light.landmarks import THAlphaHelix
from light.performance.stats import Stats

PEAK_HYDROPATHY = THAlphaHelix.PEAK_HYDROPATHY
TROUGH_HYDROPATHY = THAlphaHelix.TROUGH_HYDROPATHY


def analyzeAlphaHelix(read):
    """
    Analyze an amino acid sequence holding an alpha helix.

    @param read: An instance of C{dark.reads.AARead}.
    @return: A C{dict} with keys as follows:
        'consecutivePeaks': A C{list} of C{int} counts of the number of
            consecutive peak AAs found for each time the hydropathy graph
            goes through an overall peak. E.g., if the amino acid sequence is
            TT-PPP-TT-T-PP-T-P-PPP-T (where a 'T' is a trough hydropathy value,
            '-' is an intermediate one, and 'P' is a peak one), then
            'consecutivePeaks' will have a value of [3, 2, 4] - the count
            of the number of peak hydropathy values between troughs.
        'consecutiveTroughs': A C{list} of C{int} counts of the number of
            consecutive trough AAs found for each time the hydropathy graph
            goes through an overall trough. E.g., if the amino acid sequence is
            TT-PPP-TT-T-PP-T-P-PPP-T (where a 'T' is a trough hydropathy value,
            '-' is an intermediate one, and 'P' is a peak one), then
            'consecutiveTroughs' will have a value of [2, 3, 1, 1] - the count
            of the number of trough hydropathy values between peaks.
        'extremaCount': The C{int} number of total extreme (peak or trough)
            hydropathy values found in C{read}.
        'final': The C{int} number of amino acids at the end of C{read} that
            follow the last-found extremum, or C{None} if no extremum was
            found.
        'initial': The C{int} number of amino acids at the start of C{read}
            that before the first-found extremum, or C{None} if no extremum
            was found.
    }
    """
    # States.
    LOOKING = 0
    AWAITING_PEAK = 1
    AWAITING_TROUGH = 2

    state = LOOKING

    consecutivePeaks = []
    consecutiveTroughs = []

    for offset, aa in enumerate(read.sequence):
        try:
            hydropathy = PROPERTY_CLUSTERS[aa]['hydropathy']
        except KeyError:
            hydropathy = 0

        if state == LOOKING:
            if hydropathy >= PEAK_HYDROPATHY:
                state = AWAITING_TROUGH
                finalExtremaOffset = offset
                extremaCount = 1
                consecutivePeakCount = 1
                consecutiveTroughCount = 0
                initial = offset
            elif hydropathy <= TROUGH_HYDROPATHY:
                state = AWAITING_PEAK
                finalExtremaOffset = offset
                extremaCount = 1
                consecutivePeakCount = 0
                consecutiveTroughCount = 1
                initial = offset

        elif state == AWAITING_PEAK:
            if hydropathy >= PEAK_HYDROPATHY:
                finalExtremaOffset = offset
                state = AWAITING_TROUGH
                extremaCount += 1
                consecutiveTroughs.append(consecutiveTroughCount)
                consecutivePeakCount = 1
                consecutiveTroughCount = 0
            elif hydropathy <= TROUGH_HYDROPATHY:
                finalExtremaOffset = offset
                consecutiveTroughCount += 1
                extremaCount += 1

        elif state == AWAITING_TROUGH:
            if hydropathy <= TROUGH_HYDROPATHY:
                finalExtremaOffset = offset
                state = AWAITING_PEAK
                extremaCount += 1
                consecutivePeaks.append(consecutivePeakCount)
                consecutivePeakCount = 0
                consecutiveTroughCount = 1
            elif hydropathy >= PEAK_HYDROPATHY:
                finalExtremaOffset = offset
                consecutivePeakCount += 1
                extremaCount += 1

    if state == LOOKING:
        initial = final = None
        extremaCount = 0
        assert sum(consecutivePeaks) == 0
        assert sum(consecutiveTroughs) == 0
    else:
        final = len(read) - finalExtremaOffset - 1
        if state == AWAITING_PEAK:
            consecutiveTroughs.append(consecutiveTroughCount)
        else:
            consecutivePeaks.append(consecutivePeakCount)

    # Sanity check.
    assert extremaCount == sum(consecutivePeaks) + sum(consecutiveTroughs), (
        'Oops... %d != %r + %r (%s)' % (extremaCount, consecutivePeaks,
                                        consecutiveTroughs, read.sequence))

    return {
        'consecutivePeaks': consecutivePeaks,
        'consecutiveTroughs': consecutiveTroughs,
        'extremaCount': extremaCount,
        'final': final,
        'initial': initial,
    }


def analyzeAlphaHelices(reads):
    """
    Analyze a set of reads that contain alpha helices using
    C{analyzeAlphaHelix} and return their overall statistics.

    @return: A C{dict} with the following keys, each of which has a C{Stats}
        instance as its value:

            'length': The alpha helix lengths.
            'initial': The C{int} number of amino acids at the start of
                alpha helices that occurred before the first-found extremum.
            'final': The C{int} number of amino acids at the end of
                alpha helices that occurred after the last-found extremum.
            'extremaCount': The number of extrema found in each alpha helix.
            'consecutivePeaks': A C{list} of C{int} counts of the number of
                consecutive peak AAs found for each time the hydropathy graph
                of an alpha helix goes through an overall peak. See the
                docstring for C{analyzeAlphaHelix} above.
            'consecutiveTroughs': A C{list} of C{int} counts of the number of
                consecutive trough AAs found for each time the hydropathy graph
                of an alpha helix goes through an overall trough. See the
                docstring for C{analyzeAlphaHelix} above.
            'noExtremaLength': The lengths of the alpha helies that do not
                contain any extrema.
    """
    length = Stats('Alpha helix length')
    initial = Stats('Initial non-peak non-trough AAs')
    final = Stats('Final non-peak non-trough AAs')
    extremaCount = Stats('Number of extrema')
    consecutivePeaks = Stats('Consecutive peaks')
    consecutiveTroughs = Stats('Consecutive troughs')
    noExtremaLength = Stats('Length of helices with no extrema')

    for read in reads:
        readLen = len(read)
        length.append(readLen)
        analysis = analyzeAlphaHelix(read)

        if analysis['initial'] is None:
            # No extrema were found in the read.
            noExtremaLength.append(readLen)
            extremaCount.append(0)
        else:
            initial.append(analysis['initial'])
            final.append(analysis['final'])
            extremaCount.append(analysis['extremaCount'])
            consecutivePeaks.extend(analysis['consecutivePeaks'])
            consecutiveTroughs.extend(analysis['consecutiveTroughs'])

    return {
        'length': length,
        'initial': initial,
        'final': final,
        'extremaCount': extremaCount,
        'consecutivePeaks': consecutivePeaks,
        'consecutiveTroughs': consecutiveTroughs,
        'noExtremaLength': noExtremaLength,
    }


def selectSubstringsForAhoCorasick(alphaHelixInfo, minTruePositives=0,
                                   minTruePositiveFraction=0.0,
                                   maxSubstrings=-1):
    """
    Select alpha helix substrings for use in the Aho Corasick finder.

    @param alphaHelixInfo: An iterable of C{str}s, each having three
        space-separated fields: an alpha helix substring, the integer true
        positive count, and the integer false positive count.
    @param minTruePositives: The C{int} minimum number of true positives a
        substring must have to be considered.
    @param minTruePositiveFraction: The minimum ratio of true positives to
        overall (true + false) positives a substring must have to be
        considered.
    @param maxSubstrings: An C{int} limit on the number of substrings returned.
        If -1, no limit is applied.
    @raise ValueError: If a substring occurs more than once in
        C{alphaHelixInfo}.
    @return: A C{dict}, with C{str} keys as follows:
        fractionTooLow: The C{int} number of substrings whose true positive
            fraction was below the threshold given by
            C{minTruePositiveFraction}.
        inferior: The C{int} number of substrings that were discarded because
            they had a subsubstring whose true positive fraction was better.
        inputCount: The C{int} number of substrings given in
            C{alphaHelixInfo}.
        notEnoughTruePositives: The C{int} number of substrings whose true
            positive count was below the threshold given by
            C{minTruePositives}.
        substrings: A C{list} of (substring, (truePositives, falsePositives,
            fraction)) tuples. The tuples are sorted on true positive fraction
            (decreasing), then on substring length (increasing), and then on
            the substring itself (increasing).
    """
    substrings = {}
    inputCount = notEnoughTruePositives = fractionTooLow = inferior = 0

    for line in alphaHelixInfo:
        inputCount += 1
        (substring, truePositives, falsePositives) = line.split()
        if substring in substrings:
            raise ValueError('Substring %r found multiple times' % substring)
        truePositives = int(truePositives)
        if truePositives < minTruePositives:
            notEnoughTruePositives += 1
            continue
        falsePositives = int(falsePositives)
        fraction = truePositives / (truePositives + falsePositives)
        if fraction < minTruePositiveFraction:
            fractionTooLow += 1
            continue
        substrings[substring] = (truePositives, falsePositives, fraction)

    minLength = min(len(s) for s in substrings) if substrings else 0

    def getSubstrings(s, minLength):
        """
        Generate all substrings of a string, from longest to shortest.

        @param s: A C{str} string whose substrings are wanted.
        @param minLength: The C{int} minimal length substring to return.
        @return: A generator that yields C{str} substrings.
        """
        length = len(s)
        for substringLength in reversed(range(minLength, length)):
            numberOfSubstrings = length - substringLength + 1
            for offset in range(numberOfSubstrings):
                yield s[offset:offset + substringLength]

    # Reduce the substring set. We do this by examining substrings in turn.
    # For a given substring, if any of its substrings (i.e., the
    # sub-substrings) has an equally good or better (higher) true positive
    # fraction, we discard the substring. If not, we keep the substring.
    reducedSubstrings = []
    for substring, counts in six.iteritems(substrings):
        fraction = counts[2]
        for subsubstring in getSubstrings(substring, minLength):
            if subsubstring in substrings:
                if substrings[subsubstring][2] >= fraction:
                    # No need to keep this substring as its true positive
                    # fraction is no better than (at least) one of its
                    # substrings.
                    inferior += 1
                    break
        else:
            # No subsubstring was better. Keep this substring.
            reducedSubstrings.append((substring, counts))

    def lengthThenSubstringKey(item):
        """
        Return a tuple of a substring's length and the substring.

        @param item: A 2-tuple of
            (substring, (truePositives, falsePositives, fraction)).
        @return: A 2-C{tuple} of the substring length and the substring.
        """
        return len(item[0]), item[0]

    def fractionKey(item):
        """
        Return a substring's true positive fraction.

        @param item: A 2-tuple of
            (substring, (truePositives, falsePositives, fraction)).
        @return: The C{float} substring true positive fraction.
        """
        return item[1][2]

    reducedSubstrings.sort(key=lengthThenSubstringKey)
    reducedSubstrings.sort(key=fractionKey, reverse=True)

    if maxSubstrings > -1:
        reducedSubstrings = reducedSubstrings[:maxSubstrings]

    return {
        'fractionTooLow': fractionTooLow,
        'inferior': inferior,
        'inputCount': inputCount,
        'notEnoughTruePositives': notEnoughTruePositives,
        'substrings': reducedSubstrings,
    }
