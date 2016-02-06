from random import sample

from light.features import Landmark
from light.finder import Finder


class RandomLandmark(Finder):
    """
    A class for computing random landmarks of length 1. This landmark finder
    will generate random landmarks. The number of landarks is determined by the
    density parameter.
    """
    NAME = 'RandomLandmark'
    SYMBOL = 'RL'

    def find(self, read):
        """
        A function that returns a specified number of landmarks at random
        positions.

        @param read: An instance of C{dark.reads.AARead}.
        """
        landmarkNumber = int(len(read) * self.randomLandmarkDensity)
        offsets = sample(range(len(read)), landmarkNumber)
        for offset in offsets:
            yield Landmark(self.NAME, self.SYMBOL, offset, 1)
