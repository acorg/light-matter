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
        final = len(read) - finalExtremaOffset
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
    verbose = False

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

        if verbose:
            print('read len=%d %s id=%s seq=%r' % (
                readLen, read.id, read.sequence))

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
