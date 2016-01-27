from random import sample

from light.features import TrigPoint
from light.finder import Finder


class RandomTrigPoint(Finder):
    """
    A class for computing random trigPoints. This trigPoint finder will
    generate random trigPoints. The number of trigPoints is determined by the
    density parameter.
    """
    NAME = 'RandomTrigPoint'
    SYMBOL = 'RT'

    def find(self, read, density=None):
        """
        A function that returns a specified number of trigPoints at random
        positions.

        @param read: An instance of C{dark.reads.AARead}.
        @param density: A C{float} specifiying the density of landmarks to be
            returned.
        """
        density = density or 0.1
        landmarkNumber = int(len(read) * density)
        offsets = sample(range(len(read)), landmarkNumber)
        for offset in offsets:
            yield TrigPoint(self.NAME, self.SYMBOL, offset)
