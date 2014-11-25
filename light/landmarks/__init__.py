from light.landmarks.alpha_helix import AlphaHelix


ALL_LANDSCAPE_FINDER_CLASSES = {AlphaHelix}


def find(name):
    """
    A function to find a landmark finder.

    @param name: The name of the landmark finder to find.
    """

    for klass in ALL_LANDSCAPE_FINDER_CLASSES:
        if name == klass.__name__:
            return klass


# Default exports for 'from light.landmarks import *'
__all__ = ['find', 'ALL_LANDSCAPE_FINDER_CLASSES']
