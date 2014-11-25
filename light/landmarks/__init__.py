from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi


ALL_LANDMARK_FINDER_CLASSES = {AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi}


def find(name):
    """
    A function to find a landmark finder.

    @param name: The name of the landmark finder to find.
    """

    for klass in ALL_LANDMARK_FINDER_CLASSES:
        if name == klass.__name__:
            return klass


# Default exports for 'from light.landmarks import *'
__all__ = ['find', 'ALL_LANDMARK_FINDER_CLASSES']
