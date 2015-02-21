import re

from light.features import Landmark


class BetaTurn(object):
    """
    A class for computing statistics based on turns. Data is from:
    http://www.sciencedirect.com/science/article/pii/0022283677900948, Chou and
    Fasman (1977), Journal of Molecular Biology.
    """
    NAME = 'BetaTurn'
    SYMBOL = 'BT'
    BETA_TURN = re.compile('[NCD][PSK][NDG][WGY]')

    def find(self, read):
        """
        A function that checks if and where a turn in a sequence occurs and
        returns C{Landmark} instances.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for match in self.BETA_TURN.finditer(read.sequence):
            start = match.start()
            end = match.end()
            length = end - start
            yield Landmark(self.NAME, self.SYMBOL, start, length)
