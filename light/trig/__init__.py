from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.trig.amino_acids import AminoAcids

ALL_TRIG_FINDER_CLASSES = {Peaks, Troughs, AminoAcids}


def find(name):
    """
    A function to find a trig point finder.

    @param name: The name of the trig point finder to find.
    """

    for klass in ALL_TRIG_FINDER_CLASSES:
        if name == klass.NAME:
            return klass


# Default exports for 'from light.trig import *'
__all__ = ['find', 'ALL_TRIG_FINDER_CLASSES']
