from random import sample

from light.features import Landmark
from light.finder import Finder


class RandomLandmark(Finder):
    """
    A class for producing random landmarks of length 1. The number of landarks
    generated is determined by the database randomLandmarkDensity parameter.
    """
    NAME = 'RandomLandmark'
    SYMBOL = 'RL'

    def find(self, read):
        """
        A function that returns a specified number of landmarks at random
        positions.

        @param read: An instance of C{dark.reads.AARead}.
        """
        landmarkNumber = int(len(read) * self._dbParams.randomLandmarkDensity)
        offsets = sample(range(len(read)), landmarkNumber)
        for offset in offsets:
            yield Landmark(self.NAME, self.SYMBOL, offset, 1)
