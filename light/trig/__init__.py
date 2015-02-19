from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.trig.amino_acids import AminoAcids
from light.trig.individual_peaks import IndividualPeaks
from light.trig.individual_troughs import IndividualTroughs

ALL_TRIG_FINDER_CLASSES = {
    Peaks, Troughs, AminoAcids, IndividualPeaks, IndividualTroughs}

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


def findTrigPoints(names):
    """
    A function to find multiple trig point finders by name.

    @param names: A C{list} of C{str} name of the trig point finder classes to
        find, or C{None} (equivalent to an empty list of names).
    @raise ValueError: If any name cannot be found.
    @return: A C{list} of classes that were successfully found.
    """
    found = []
    unknown = []
    if names:
        for name in names:
            klass = findTrigPoint(name)
            if klass:
                found.append(klass)
            else:
                unknown.append(name)
    if unknown:
        raise ValueError('Unknown trig point finder%s: %s.' % (
            's' if len(unknown) > 1 else '', ', '.join(unknown)))
    else:
        return found


# Default exports for 'from light.trig import *'
__all__ = ['findTrigPoint', 'findTrigPoints', 'ALL_TRIG_FINDER_CLASSES',
           'DEFAULT_TRIG_FINDER_CLASSES']
