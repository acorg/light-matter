from unittest import TestCase
from os.path import dirname, join
from json import loads

from dark.reads import AARead

import light
from light.features import Landmark
from light.landmarks.cluster_alpha_helix import ClusterAlphaHelix


class TestClusterAlphaHelix(TestCase):
    """
    Tests for the Landmark.ClusterAlphaHelix class.
    """

    def testFindWithoutCluster(self):
        """
        The find method must return an empty generator when no alpha helix
        cluster is present.
        """
        read = AARead('id', 'FFFFFFFFFFFFFFFFFFFFFFFFFFF')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFindOneClusterBeginning(self):
        """
        The find method must find one cluster at the beginning of a sequence.
        """
        read = AARead('id', 'KKAHFFFFFFFFF')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131')],
                         result)

    def testFindOneClusterMiddle(self):
        """
        The find method must find a cluster in the middle of a sequence.
        """
        read = AARead('id', 'AAAAAAAAKKAHAAAAAAAA')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 8, 4, '1131')],
                         result)

    def testFindOneClusterEnd(self):
        """
        The find method must find a cluster at the end of a sequence.
        """
        read = AARead('id', 'AAAAAAAAKKAH')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 8, 4, '1131')],
                         result)

    def testFindTwoIdenticalAAClustersNotAdjacent(self):
        """
        The find method must find repeating identical clusters that come from
        the same AA sequence.
        """
        read = AARead('id', 'KKAHAAAAAAAAKKAH')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 12, 4, '1131')],
                         result)

    def testFindTwoIdenticalClustersNonIdenticalAANotAdjacent(self):
        """
        The find method must find repeating identical clusters that come from
        a different AA sequence.
        """
        read = AARead('id', 'KKAHAAAAAAAARKMH')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 12, 4, '1131')],
                         result)

    def testFindTwoIdenticalAAClustersAdjacent(self):
        """
        The find method must find two identical clusters with the same AA
        sequence.
        """
        read = AARead('id', 'KKAHFFFKKAH')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 7, 4, '1131')],
                         result)

    def testFindTwoIdenticalAAClustersNonIdenticalAAAdjacent(self):
        """
        The find method must find two identical clusters with a different AA
        sequence.
        """
        read = AARead('id', 'KKAHRKMH')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 4, 4, '1131')],
                         result)

    def testThreeDifferentClustersTwoOccurences(self):
        """
        The find method must find three different clusters with two occurences
        each.
        """
        read = AARead('id', 'KKAHGHLVKVIGKKAHGHLVKVIG')
        finder = ClusterAlphaHelix()
        result = list(finder.find(read))
        self.assertEqual([Landmark('ClusterAlphaHelix', 'CAH', 0, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 4, 4, '2144'),
                          Landmark('ClusterAlphaHelix', 'CAH', 8, 4, '1442'),
                          Landmark('ClusterAlphaHelix', 'CAH', 12, 4, '1131'),
                          Landmark('ClusterAlphaHelix', 'CAH', 16, 4, '2144'),
                          Landmark('ClusterAlphaHelix', 'CAH', 20, 4, '1442')],
                         result)

    def testClustersAreUnique(self):
        """
        All clusters in the cluster database must be unique.
        """
        filename = join(dirname(light.__file__), '..', 'data',
                        'cluster_alpha_helix_3.25.json')
        database = []
        append = database.append
        with open(filename) as fp:
            for line in fp:
                clusterLine = loads(line)
                append(clusterLine['cluster'])
        # The 17 is based on the cluster database of March 7, 2016.
        self.assertEqual(17, len(database))
        self.assertEqual(len(database), len(set(database)))
