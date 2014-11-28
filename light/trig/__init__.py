from light.trig.peaks import Peaks
from light.trig.troughs import Troughs

ALL_TRIG_FINDER_CLASSES = {Peaks, Troughs}


def find(name):
    """
    A function to find a landmark finder.

    @param name: The name of the landmark finder to find.
    """

    for klass in ALL_TRIG_FINDER_CLASSES:
        if name == klass.NAME:
            return klass


# Default exports for 'from light.landmarks import *'
__all__ = ['find', 'ALL_TRIG_FINDER_CLASSES']
