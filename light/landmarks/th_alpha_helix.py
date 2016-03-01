from dark.aa import PROPERTY_CLUSTERS

from light.distance import scaleLog
from light.features import Landmark
from light.finder import Finder


class THAlphaHelix(Finder):
    """
    A shoot-first-ask-questions-later trigger-happy alpha helix finder.
    """
    NAME = 'THAlphaHelix'
    SYMBOL = 'THA'

    # The cut-off values below are all inclusive. We may want to move
    # some/all of these into light.parameters.DatabaseParameters for access
    # via self._dbParams, but for now they are constants.

    # Maximum distances (i.e., number of AAs) that we allow when waiting to
    # encounter a hydropathy peak or trough.
    MAX_WAIT_FOR_PEAK = 6
    MAX_WAIT_FOR_TROUGH = 6

    # Once we see a peak (trough) AA, we then wait to see a trough (peak)
    # to continue the helix. But in the meantime we may see another peak
    # (trough). The following two place a limit on how many consecutive
    # peaks (troughs) we tolerate before giving up on trying to extend the
    # helix. The statistics after each give the values found via analyzing
    # all PDB alpha helices.
    MAX_CONSECUTIVE_PEAKS = 5    # Mean: 1.96, Median: 2, SD: 1.42
    MAX_CONSECUTIVE_TROUGHS = 4  # Mean: 1.61, Median: 1, SD: 0.99

    # The AA (cluster) hydropathy values considered large/small enough to
    # be a peak/trough.
    PEAK_HYDROPATHY = 3
    TROUGH_HYDROPATHY = 1

    MIN_HELIX_LENGTH = 7

    # The minimum total number of peak/trough AAs that must be present in
    # an alpha helix.
    MIN_EXTREMA_COUNT = 3

    def find(self, read):
        """
        Find possible alpha helices in a sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """

        # States.
        LOOKING = 0
        AWAITING_PEAK = 1
        AWAITING_TROUGH = 2

        state = LOOKING

        for offset, aa in enumerate(read.sequence):
            try:
                hydropathy = PROPERTY_CLUSTERS[aa]['hydropathy']
            except KeyError:
                hydropathy = 0

            if state == LOOKING:
                if hydropathy >= self.PEAK_HYDROPATHY:
                    state = AWAITING_TROUGH
                    startOffset = lastGoodOffset = offset
                    maxWait = self.MAX_WAIT_FOR_TROUGH
                    extremaCount = 1
                    peakCount = 1
                elif hydropathy <= self.TROUGH_HYDROPATHY:
                    state = AWAITING_PEAK
                    startOffset = lastGoodOffset = offset
                    maxWait = self.MAX_WAIT_FOR_PEAK
                    extremaCount = 1
                    troughCount = 1

            elif offset - lastGoodOffset >= maxWait:
                # We've waited too long to extend the current helix.
                helix = self._isHelix(startOffset, lastGoodOffset,
                                      extremaCount)
                if helix:
                    yield helix
                state = LOOKING

            elif state == AWAITING_PEAK:
                if hydropathy >= self.PEAK_HYDROPATHY:
                    # Found the peak we were waiting for.
                    state = AWAITING_TROUGH
                    lastGoodOffset = offset
                    maxWait = self.MAX_WAIT_FOR_TROUGH
                    extremaCount += 1
                    peakCount = 1
                elif hydropathy <= self.TROUGH_HYDROPATHY:
                    # We found an additional (consecutive) trough.
                    if troughCount == self.MAX_CONSECUTIVE_TROUGHS:
                        helix = self._isHelix(startOffset, lastGoodOffset,
                                              extremaCount)
                        if helix:
                            yield helix
                        # There is a subtlety here. The switch back to
                        # LOOKING state means the trough we just found is
                        # not treated as the potential start of a new alpha
                        # helix. The alternative is to go straight into
                        # AWAITING_PEAK state (and set variables as in the
                        # 'if state == LOOKING' code above).
                        state = LOOKING
                    else:
                        troughCount += 1
                        lastGoodOffset = offset

            elif state == AWAITING_TROUGH:
                if hydropathy <= self.TROUGH_HYDROPATHY:
                    # Found the trough we were waiting for.
                    state = AWAITING_PEAK
                    lastGoodOffset = offset
                    maxWait = self.MAX_WAIT_FOR_PEAK
                    extremaCount += 1
                    troughCount = 1
                elif hydropathy >= self.PEAK_HYDROPATHY:
                    # We found an additional (consecutive) peak.
                    if peakCount == self.MAX_CONSECUTIVE_PEAKS:
                        helix = self._isHelix(startOffset, lastGoodOffset,
                                              extremaCount)
                        if helix:
                            yield helix
                        # There is a subtlety here. The switch back to
                        # LOOKING state means the peak we just found is not
                        # treated as the potential start of a new alpha
                        # helix. The alternative is to go straight into
                        # AWAITING_TROUGH state (and set variables as in
                        # the 'if state == LOOKING' code above).
                        state = LOOKING
                    else:
                        peakCount += 1
                        lastGoodOffset = offset

        # Yield the final helix, if any.
        if state != LOOKING:
            helix = self._isHelix(startOffset, lastGoodOffset, extremaCount)
            if helix:
                yield helix

    def _isHelix(self, startOffset, lastGoodOffset, extremaCount):
        """
        Helper function to check whether the current subsequence (beginning at
        C{startOffset} in a read can be emitted by C{find} as a possible helix.

        @param startOffset: The C{int} offset where the possible helix started
            in the read.
        @param lastGoodOffset: The C{int} offset where we encountered the last
            AA (an extremum) suspected of being in the helix.
        @param extremaCount: The C{int} number of extrema found in this
            possible helix.
        @return: A C{Landmark} instance if the subsequence is of sufficient
            length and has enough extrema to qualify as a helix. Else C{None}.
        """
        helixLength = lastGoodOffset - startOffset + 1
        if (helixLength >= self.MIN_HELIX_LENGTH and
                extremaCount >= self.MIN_EXTREMA_COUNT):
            return Landmark(
                self.NAME, self.SYMBOL, startOffset, helixLength,
                scaleLog(helixLength, self._dbParams.featureLengthBase))
