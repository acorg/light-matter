from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.individual_peaks import IndividualPeaks


class TestIndividualPeaks(TestCase):
    """
    Tests for the TrigPoint.IndividualPeaks class.
    """
    def testFindWithoutIndividualPeak(self):
        """
        The find method must return an empty generator when no individual peak
        is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        peaks = IndividualPeaks()
        result = list(peaks.find(read))
        self.assertEqual([], result)

    def testTwoIndividualPeaks(self):
        """
        The find method must find two individual peaks.
        """
        read = AARead('id', 'ATAAAAATAAA')
        peaks = IndividualPeaks()
        result = list(peaks.find(read))
        self.assertEqual([TrigPoint('IndividualPeaks', 'I', 1),
                          TrigPoint('IndividualPeaks', 'I', 7)], result)

    def testNoIndividualPeakAtBeginning(self):
        """
        The find method must not find an individual peak at the beginning of
        the sequence.
        """
        read = AARead('id', 'TAAAAAAAAA')
        peaks = IndividualPeaks()
        result = list(peaks.find(read))
        self.assertEqual([], result)

    def testNoIndividualPeakAtEnd(self):
        """
        The find method must not find an individual peak at the end of the
        sequence.
        """
        read = AARead('id', 'AAAAAAT')
        peaks = IndividualPeaks()
        result = list(peaks.find(read))
        self.assertEqual([], result)
