from light.landmarks.alpha_helix import AlphaHelix


def find(name):
    """
    A function to find a landmark finder.

    @param name: The name of the landmark finder to find.
    """

    if name == AlphaHelix.__name__:
        return AlphaHelix
    else:
        return None


# Default exports for 'from light.landmarks import *'
__all__ = ['find']
