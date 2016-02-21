from __future__ import absolute_import

from random import sample

from light.features import TrigPoint
from light.finder import Finder


class RandomTrigPoint(Finder):
    """
    A class for computing random trigPoints. The number of trig points
    produced is determined by the database randomTrigPointDensity parameter.
    """
    NAME = 'RandomTrigPoint'
    SYMBOL = 'RT'

    def find(self, read):
        """
        A function that returns a specified number of trigPoints at random
        positions.

        @param read: An instance of C{dark.reads.AARead}.
        """
        landmarkNumber = int(len(read) * self._dbParams.randomTrigPointDensity)
        offsets = sample(range(len(read)), landmarkNumber)
        for offset in offsets:
            yield TrigPoint(self.NAME, self.SYMBOL, offset)
