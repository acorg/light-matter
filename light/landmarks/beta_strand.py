import re

from light.features import Landmark


class BetaStrand(object):
    """
    Simplistic and permissive (i.e., this *will* produce false positives)
    identification of amino acid sequences that hopefully correspond to Beta
    strands.

    Note that false positives may not be a problem. Even if a long sequence of
    AAs with high beta strand propensities is not in fact a beta strand, their
    existence may still be statistical relevant for matching purposes.

    This is based on simple beta strand amino acid propensities, as described
    in "Dependence of a-helical and b-sheet amino acid propensities on the
    overall protein fold type" (Kazuo Fujiwara, Hiromi Toda and Masamichi
    Ikeguchi, 2012) at http://www.biomedcentral.com/1472-6807/12/18

    The AAs permitted to occur in what we consider a beta strand (see
    BETASTRAND, below) are quite permissive.  See Table 2 of the paper
    mentioned above. If an AA has *either* a buried or exposed propensity of
    1.0 (or close) it is included.

    The minimal allowed length of a beta strand is 3. Taken from Figure 1 of
    "Position-specific propensities of amino acids in the b-strand" (Nicholus
    Bhattacharjee, Parbati Biswas, 2010) at
    http://www.biomedcentral.com/1472-6807/10/29  We do not impose an upper
    limit on length (despite no known beta strands being over 15 AAs) as this
    would rule out beta strands that occur within a longer sequence of
    acceptable AAs.
    """
    NAME = 'BetaStrand'
    SYMBOL = 'S'
    BETA_STRAND_AAs = 'VILMCFYWQTHKR'
    MIN_LENGTH = 3
    _BETA_STRAND_RE = re.compile('[%s]{%d,}' % (BETA_STRAND_AAs, MIN_LENGTH))

    def find(self, read):
        """
        Find possible beta strands in a sequence (see above comments re
        approach and false positives, etc).

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for match in self._BETA_STRAND_RE.finditer(read.sequence):
            start = match.start()
            end = match.end()
            length = end - start
            yield Landmark(self.NAME, self.SYMBOL, start, length, length)
