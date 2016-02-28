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

    # Maximum distances (i.e., number of AAs) that we allow when waiting to
    # encounter a hydropathy peak or trough.
    MAX_WAIT_FOR_PEAK = 6
    MAX_WAIT_FOR_TROUGH = 6

    # Cut-off values. These are inclusive (i.e., they are all tested
    # using <= or >=).
    PEAK_HYDROPATHY = 3
    TROUGH_HYDROPATHY = 1
    MIN_HELIX_LENGTH = 7
    MIN_EXTREMA_COUNT = 3

    def find(self, read):
        """
        Find potential alpha helices in a sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        featureLengthBase = self._dbParams.featureLengthBase

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
                elif hydropathy <= self.TROUGH_HYDROPATHY:
                    state = AWAITING_PEAK
                    startOffset = lastGoodOffset = offset
                    maxWait = self.MAX_WAIT_FOR_PEAK
                    extremaCount = 1

            elif offset - lastGoodOffset >= maxWait:
                # We've waited too long. Yield the helix up to the last
                # good offset, if it was long enough and enough
                # peaks/troughs were seen.
                helixLength = lastGoodOffset - startOffset + 1
                if (helixLength >= self.MIN_HELIX_LENGTH and
                        extremaCount >= self.MIN_EXTREMA_COUNT):
                    yield Landmark(
                        self.NAME, self.SYMBOL, startOffset, helixLength,
                        scaleLog(helixLength, featureLengthBase))
                state = LOOKING

            elif state == AWAITING_PEAK:
                if hydropathy >= self.PEAK_HYDROPATHY:
                    lastGoodOffset = offset
                    state = AWAITING_TROUGH
                    maxWait = self.MAX_WAIT_FOR_TROUGH
                    extremaCount += 1
                elif hydropathy <= self.TROUGH_HYDROPATHY:
                    # We found another trough. Let's remember it.
                    lastGoodOffset = offset

            elif state == AWAITING_TROUGH:
                if hydropathy <= self.TROUGH_HYDROPATHY:
                    lastGoodOffset = offset
                    state = AWAITING_PEAK
                    maxWait = self.MAX_WAIT_FOR_PEAK
                    extremaCount += 1
                elif hydropathy >= self.PEAK_HYDROPATHY:
                    # We found another peak. Let's remember it.
                    lastGoodOffset = offset

        # Yield the final helix, if any.
        if state != LOOKING:
            helixLength = lastGoodOffset - startOffset + 1
            if (helixLength >= self.MIN_HELIX_LENGTH and
                    extremaCount >= self.MIN_EXTREMA_COUNT):
                yield Landmark(
                    self.NAME, self.SYMBOL, startOffset, helixLength,
                    scaleLog(helixLength, featureLengthBase))
