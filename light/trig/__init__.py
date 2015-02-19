from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.trig.amino_acids import AminoAcids
from light.trig.polarity_peaks import PolarityPeaks
from light.trig.individual_peaks import IndividualPeaks
from light.trig.individual_troughs import IndividualTroughs


ALL_TRIG_FINDER_CLASSES = {Peaks, Troughs, AminoAcids, IndividualPeaks,
                           IndividualTroughs, PolarityPeaks}

DEFAULT_TRIG_FINDER_CLASSES = {
    Peaks, Troughs, AminoAcids}


def findTrigPoint(name):
    """
    A function to find a trig point finder.

    @param name: The C{str} name of the trig point finder class to find.
    @return: The found class, or C{None}.
    """

    for klass in ALL_TRIG_FINDER_CLASSES:
        if name == klass.NAME:
            return klass


# Default exports for 'from light.trig import *'
__all__ = ['findTrigPoint', 'ALL_TRIG_FINDER_CLASSES',
           'DEFAULT_TRIG_FINDER_CLASSES']
