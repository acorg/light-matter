import re
from json import loads
from os.path import dirname, join

from dark.aa import PROPERTY_CLUSTERS

import light
from light.features import Landmark
from light.finder import Finder


def _loadDatabase():
    """
    Read and convert the JSON file of alpha helix cluster patterns.

    The set of alpha helix cluster patterns was constructed as follows:
    - All 'H' alpha helices were extracted from
      http://www.rcsb.org/pdb/files/ss.txt.
    - The alpha helices were converted to cluster patterns using the clusters
      in C{dark.aa}. Duplicates and clusters containing 0 were removed.
    - True and false positives were evaluated using all sequences in
      http://www.rcsb.org/pdb/files/ss.txt.
    - All patterns that had a true positive to false positive rate smaller or
      equal of 3.25 were included in this finder.

    @return: A C{list} of C{dict}s, each containing 1) a 'name' key
        giving the C{str} alpha helix cluster and 2) a 'regex' key with
        a compiled regex value.
    """
    filename = join(dirname(light.__file__), '..', 'data',
                    'cluster_alpha_helix_3.25.json')
    database = []
    append = database.append
    with open(filename) as fp:
        for line in fp:
            clusterLine = loads(line)
            regex = re.compile(clusterLine['cluster'])
            append({
                'cluster': clusterLine['cluster'],
                'regex': regex,
            })
    return database


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
    clusterRead = ''
    for aa in read:
        try:
            clusterRead += str(PROPERTY_CLUSTERS[aa]['hydropathy'])
        except KeyError:
            clusterRead += str(0)
    return clusterRead


class ClusterAlphaHelix(Finder):
    """
    A class for finding prosite motifs.
    """
    NAME = 'ClusterAlphaHelix'
    SYMBOL = 'CAH'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) cluster
        alpha helix from a sequence

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        clusterRead = translateToCluster(read.sequence)
        for cluster in _DATABASE:
            for match in cluster['regex'].finditer(clusterRead):
                start = match.start()
                length = match.end() - start
                yield Landmark(self.NAME, self.SYMBOL, start, length,
                               cluster['cluster'])
