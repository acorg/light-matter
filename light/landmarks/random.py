from random import randint

from light.features import Landmark, Finder


class RandomLandmark(Finder):
    """
    A class for computing random landmarks of length 1. This landmark finder
    will generate random landmarks. The number of landarks is determined by the
    density parameter.
    """
    NAME = 'RandomLandmark'
    SYMBOL = 'RL'

    def find(self, read, density=None):
        """
        A function that returns a specified number of landmarks at random
        positions.

        @param read: An instance of C{dark.reads.AARead}.
        @param density: A C{float} specifiying the density of landmarks to be
            returned.
        """
        density = density or 0.1
        landmarkNumber = int(len(read) * density)
        offsets = [randint(0, len(read)) for i in range(landmarkNumber)]
        for offset in offsets:
            yield Landmark(self.NAME, self.SYMBOL, offset, 1)
