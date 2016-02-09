from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.trig.amino_acids import AminoAcids
from light.trig.individual_peaks import IndividualPeaks
from light.trig.individual_troughs import IndividualTroughs
from light.trig.random import RandomTrigPoint
from light.trig.volume import Volume

ALL_TRIG_CLASSES = [
    AminoAcids, IndividualPeaks, IndividualTroughs, Peaks, Troughs, Volume]

ALL_TRIG_CLASSES_EVEN_BAD_ONES = ALL_TRIG_CLASSES + [RandomTrigPoint]

DEFAULT_TRIG_CLASSES = [AminoAcids]

_HASHKEY_TO_NAME = dict((cls.SYMBOL, cls.NAME) for cls in ALL_TRIG_CLASSES)


def findTrigPoint(name):
    """
    A function to find a trig point finder.

    @param name: The C{str} name of the trig point finder class to find.
    @return: The found class, or C{None}.
    """

    for klass in ALL_TRIG_CLASSES:
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


def trigNameFromHashkey(hashkey):
    """
    Get the name of a trig point class given part of a hashkey that was
    created for it (in features.TrigPoint.hashkey).

    @param hashkey: A C{str} hashkey for the class.
    @return: A C{str} trig point class name, or C{None} if the hashkey
        cannot be found.
    """
    return _HASHKEY_TO_NAME.get(hashkey)


# Default exports for 'from light.trig import *'
__all__ = ['findTrigPoint', 'findTrigPoints', 'trigNameFromHashkey',
           'ALL_TRIG_CLASSES', 'DEFAULT_TRIG_CLASSES']
