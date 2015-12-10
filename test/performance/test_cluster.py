import six
from unittest import TestCase
import numpy as np
from collections import Counter
from sklearn.cluster import AffinityPropagation, KMeans

from dark.reads import Reads, AARead

from light.performance.cluster import (
    ClusterAnalysis, AffinityPropagationAnalysis, KMeansAnalysis)


class TestClusterAnalysis(TestCase):
    """
    Tests for the light.performance.cluster.ClusterAnalysis class.
    """

    def testNoReads(self):
        """
        The ClusterAnalysis initializer must raise ValueError if no reads
        are provided.
        """
        reads = Reads()
        error = 'No sequences were found in'
        six.assertRaisesRegex(self, ValueError, error, ClusterAnalysis, reads,
                              {}, landmarkNames=['AlphaHelix'])

    def testMissingLabelNoDefault(self):
        """
        The ClusterAnalysis initializer must raise ValueError if a read is
        given without a corresponding label and there is no default label.
        """
        reads = Reads()
        reads.add(AARead('id', 'MMM'))
        error = "Read 'id' has no corresponding label"
        six.assertRaisesRegex(self, ValueError, error, ClusterAnalysis, reads,
                              {}, landmarkNames=['AlphaHelix'])

    def testDefaultLabel(self):
        """
        The ClusterAnalysis initializer must assign the default label to any
        read without a corresponding label.
        """
        reads = Reads()
        reads.add(AARead('id', 'MMM'))
        ca = ClusterAnalysis(reads, {}, defaultLabel=5,
                             landmarkNames=['AlphaHelix'])
        self.assertEqual([5], ca.trueLabels)

    def testOneReadAffinity(self):
        """
        The ClusterAnalysis initializer must assign a 1.0 affinity to a single
        read.
        """
        reads = Reads()
        reads.add(AARead('id', 'MMM'))
        ca = ClusterAnalysis(reads, {}, defaultLabel=5,
                             landmarkNames=['AlphaHelix'])
        self.assertEqual([[1.0]], ca.affinity)

    def testTwoReadsAffinity(self):
        """
        The ClusterAnalysis initializer must assign a 1.0 affinity to a read
        with itself and a 0.0 affinity to a second read that the first has
        nothing in common with.
        """
        reads = Reads()
        reads.add(AARead('id1', 'MMM'))
        reads.add(AARead('id2', 'MMM'))
        ca = ClusterAnalysis(reads, {}, defaultLabel=5,
                             landmarkNames=['AlphaHelix'])
        expected = np.ones((2, 2))
        expected[0, 1] = expected[1, 0] = 0.0
        for row in range(2):
            for col in range(2):
                self.assertEqual(expected[row, col], ca.affinity[row, col])

    def testAffinityFromOffsetDeltasEmpty(self):
        """
        The ClusterAnalysis.affinityFromOffsetDeltas function must return
        0.0 when two empty offset counters are passed.
        """
        result = ClusterAnalysis.affinityFromOffsetDeltas(
            Counter(), Counter())
        self.assertEqual(0.0, result)

    def testAffinityFromOffsetDeltasNothingInCommon(self):
        """
        The ClusterAnalysis.affinityFromOffsetDeltas function must return
        0.0 when two sets of offset deltas have nothing in common.
        """
        result = ClusterAnalysis.affinityFromOffsetDeltas(
            Counter({46: 2, 47: 3}),
            Counter({33: 4, 38: 1}))
        self.assertEqual(0.0, result)

    def testAffinityFromOffsetDeltasIdentical(self):
        """
        The ClusterAnalysis.affinityFromOffsetDeltas function must return
        1.0 when two sets of offset deltas are identical.
        """
        deltas = Counter({46: 2, 47: 3})
        result = ClusterAnalysis.affinityFromOffsetDeltas(deltas, deltas)
        self.assertEqual(1.0, result)

    def testAffinityFromOffsetDeltasOneOfThreeOverlaps(self):
        """
        The ClusterAnalysis.affinityFromOffsetDeltas function must return 1/3
        when two sets of offset deltas have one thing in common and each
        has one extra delta not in common.
        """
        result = ClusterAnalysis.affinityFromOffsetDeltas(
            Counter({46: 1, 47: 1}),
            Counter({46: 1, 48: 1}))
        self.assertEqual(1.0 / 3.0, result)

    def testAffinityFromOffsetDeltasThreeOfFiveOverlaps(self):
        """
        The ClusterAnalysis.affinityFromOffsetDeltas function must return 0.6
        when two sets of offset deltas have three things in common and each
        has one extra delta not in common.
        """
        result = ClusterAnalysis.affinityFromOffsetDeltas(
            Counter({46: 2, 47: 1, 48: 1}),
            Counter({46: 2, 47: 1, 49: 1}))
        self.assertEqual(0.6, result)


class TestAffinityPropagationAnalysis(TestCase):
    """
    The the AffinityPropagationAnalysis class.
    """

    def testAffinityPropagationAttribute(self):
        """
        The AffinityPropagationAnalysis class must assign an
        AffinityPropagation instance to itself after clustering.
        """
        reads = Reads()
        reads.add(AARead('id1', 'MMM'))
        reads.add(AARead('id2', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertTrue(isinstance(ap.affinityPropagation,
                                   AffinityPropagation))

    def testTwoReadsNothingInCommonTwoClusters(self):
        """
        The AffinityPropagationAnalysis class must find two clusters when given
        two reads that have nothing in common.
        """
        reads = Reads()
        reads.add(AARead('id1', 'MMM'))
        reads.add(AARead('id2', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertEqual(2, ap.nClusters)

    def testTwoReadsNothingInCommonLabels(self):
        """
        The AffinityPropagationAnalysis class must find two clusters when given
        two reads that have nothing in common, and assign them the expected
        labels.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRF'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertEqual([0, 1], list(ap.clusterLabels))

    def testThreeReadsTwoIdenticalTwoClusters(self):
        """
        The AffinityPropagationAnalysis class must find two clusters when given
        two reads that have all deltas in common and one that has no deltas.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertEqual(2, ap.nClusters)

    def testThreeReadsTwoIdenticalLabels(self):
        """
        The AffinityPropagationAnalysis class must find two clusters when given
        two reads that have all deltas in common and one that has no deltas,
        and assign them the expected labels.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertEqual([0, 0, 1], list(ap.clusterLabels))

    def testFourReadsTwoIdenticalPrintWithoutClustering(self):
        """
        The AffinityPropagationAnalysis class must print the expected result
        when given two reads that have all deltas in common and two that have
        no deltas, even if its cluster method has not been called.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        reads.add(AARead('id4', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        # Print without calling ap.cluster()
        self.assertEqual(
            'Estimated number of clusters: 2\n'
            'Homogeneity: 1.000\n'
            'Completeness: 0.000\n'
            'V-measure: 0.000\n'
            'Adjusted Rand Index: 0.000\n'
            'Adjusted Mutual Information: 0.000\n'
            'Silhouette Coefficient: 0.667',
            ap.print_())

    def testFourReadsTwoIdenticalPrint(self):
        """
        The AffinityPropagationAnalysis class must print the expected result
        when given two reads that have all deltas in common and two that have
        no deltas.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        reads.add(AARead('id4', 'MMM'))
        ap = AffinityPropagationAnalysis(reads, {}, defaultLabel=0,
                                         landmarkNames=['AlphaHelix'])
        ap.cluster()
        self.assertEqual(
            'Estimated number of clusters: 2\n'
            'Homogeneity: 1.000\n'
            'Completeness: 0.000\n'
            'V-measure: 0.000\n'
            'Adjusted Rand Index: 0.000\n'
            'Adjusted Mutual Information: 0.000\n'
            'Silhouette Coefficient: 0.667',
            ap.print_())


class TestKMeansAnalysis(TestCase):
    """
    The the KMeansAnalysis class.
    """

    def testKMeansAttribute(self):
        """
        The KMeansAnalysis class must assign an KMeans instance to itself
        after clustering.
        """
        reads = Reads()
        reads.add(AARead('id1', 'MMM'))
        reads.add(AARead('id2', 'MMM'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(1)
        self.assertTrue(isinstance(km.kMeans, KMeans))

    def testOneReadClusterCount(self):
        """
        The KMeansAnalysis class must find one cluster when given one read,
        and set the value of nClusters to 1.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRF'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(1)
        self.assertEqual(1, km.nClusters)

    def testTwoReadsNothingInCommonClusterCount(self):
        """
        The KMeansAnalysis class must find two clusters when given two reads
        that have nothing in common, and set the value of nClusters to 2.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRF'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(2)
        self.assertEqual(2, km.nClusters)

    def testTwoReadsNothingInCommonLabels(self):
        """
        The KMeansAnalysis class must find two clusters when given two reads
        that have nothing in common, and assign them differing labels.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRF'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(2)
        result = list(km.clusterLabels)
        self.assertNotEqual(result[0], result[1])

    def testThreeReadsTwoIdenticalLabels(self):
        """
        The KMeansAnalysis class must find two clusters when given two reads
        that have all deltas in common and one that has no deltas, and assign
        the the first two the same label and the third a different label.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(2)
        result = list(km.clusterLabels)
        self.assertEqual(result[0], result[1])
        self.assertNotEqual(result[0], result[2])

    def testFourReadsTwoIdenticalLabels(self):
        """
        The KMeansAnalysis class must find two clusters when given two reads
        that have all deltas in common and two that have no deltas, and assign
        the the first two the same label, the third and the fourth the same
        label, and the labels of the two clusters must differ.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        reads.add(AARead('id3', 'MMM'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(2)
        result = list(km.clusterLabels)
        self.assertEqual(result[0], result[1])
        self.assertEqual(result[2], result[3])
        self.assertNotEqual(result[0], result[3])

    def testFourReadsTwoIdenticalPrintWithoutClustering(self):
        """
        The KMeansAnalysis class must raise RuntimeError if print_ is called
        but its cluster method has not been called.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        reads.add(AARead('id4', 'MMM'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        six.assertRaisesRegex(self, RuntimeError, '', km.print_)

    def testFourReadsTwoIdenticalPrint(self):
        """
        The KMeansAnalysis class must print the expected result
        when given two reads that have all deltas in common and two that have
        no deltas.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id3', 'MMM'))
        reads.add(AARead('id4', 'MMM'))
        km = KMeansAnalysis(reads, {}, defaultLabel=0,
                            landmarkNames=['AlphaHelix'])
        km.cluster(2)
        self.assertEqual(
            'Homogeneity: 1.000\n'
            'Completeness: 0.000\n'
            'V-measure: 0.000\n'
            'Adjusted Rand Index: 0.000\n'
            'Adjusted Mutual Information: 0.000\n'
            'Silhouette Coefficient: 0.667',
            km.print_())
