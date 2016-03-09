import re
from json import loads
from os.path import dirname, join

from dark.aa import PROPERTY_CLUSTERS

import light
from light.features import Landmark
from light.finder import Finder


# The set of hydropathy cluster alpha helix patterns was constructed as
# follows:
# - All 'H' alpha helices were extracted from
#   http://www.rcsb.org/pdb/files/ss.txt.
# - The alpha helices were converted to cluster patterns using the
#   hydrophobicity clusters in C{dark.aa}. Duplicates and clusters containing
#   0 were removed.
# - True and false positives were evaluated for a small subset of hydropathy
#   cluster alpha helices using all sequences in
#   http://www.rcsb.org/pdb/files/ss.txt.
# - All patterns with a length between 4 - 6 that had a true positive to false
#   positive rate less than or equal to 3.25 were included in this finder.


def _loadDatabase():
    """
    Read and convert the JSON file of alpha helix cluster patterns.

    @return: A compiled regular expression of hydropathy cluster alpha helices.
    """
    filename = join(dirname(light.__file__), '..', 'data',
                    'cluster_alpha_helix_3.25.json')
    database = []
    append = database.append
    with open(filename) as fp:
        for line in fp:
            clusterLine = loads(line)
            append(clusterLine['cluster'])
    regexDatabase = re.compile('|'.join(database))
    return regexDatabase


# Read the database (just once) into a singleton for use in the
# ClusterAlphaHelix.find method.
_DATABASE = _loadDatabase()


def translateToCluster(read):
    """
    Convert a C{str} amino acid sequence into a string of hydropathy property
    clusters.

    @param read: A C{str} amino acid sequence.

    @return: A C{str} of hydropathy property clusters for the input amino acid
        sequence.
    """
    clusterRead = []
    append = clusterRead.append

    for aa in read:
        try:
            append(PROPERTY_CLUSTERS[aa]['hydropathy'])
        except KeyError:
            append(0)
    return ''.join(map(str, clusterRead))


class ClusterAlphaHelix(Finder):
    """
    A class for finding alpha helices based on hydropathy clusters.
    """
    NAME = 'ClusterAlphaHelix'
    SYMBOL = 'CAH'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) alpha
        helices based on hydropathy clusters from a sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        clusterRead = translateToCluster(read.sequence)
        for match in _DATABASE.finditer(clusterRead):
            cluster = match.group()
            start = match.start()
            length = match.end() - start
            yield Landmark(self.NAME, self.SYMBOL, start, length, cluster)
