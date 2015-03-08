import sys
from collections import Counter
import numpy as np

from sklearn.cluster import AffinityPropagation, KMeans
from sklearn.metrics import (
    homogeneity_score, completeness_score, v_measure_score,
    adjusted_rand_score, adjusted_mutual_info_score, silhouette_score)

from dark.reads import AARead
from dark.fasta import FastaReads

from light.distance import scale
from light.database import DatabaseSpecifier


class ClusterAnalysis(object):

    def __init__(self, fastaFile, labels, defaultLabel=None, **kwargs):
        """
        Base class for using cluster analysis to evaluate how well various
        feature finders and database parameter settings can separate a set of
        sequences. The clustering is based on feature offset deltas.

        @param fastaFile: Either A C{str} filename of sequences to consider or
            a C{light.reads.Reads} instance.
        @param labels: A C{dict} with a label for each sequence id in
            fastaFile. These are the known categories of each sequence.
        @param defaultLabel: If not C{None}, a label to use for reads whose ids
            are not present in C{labels}. If C{None} and a read id has no label
            a ValueError is raised.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        @raises ValueError: If the id of a read is not in labels and no default
            label has been set, or if there are no reads in C{fastaFile}.
        """
        if isinstance(fastaFile, basestring):
            reads = FastaReads(fastaFile, readClass=AARead)
        else:
            reads = fastaFile
        database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        allOffsetDeltas = []
        trueLabels = []

        for read in reads:
            trueLabel = labels.get(read.id, defaultLabel)
            if trueLabel is None:
                raise ValueError('Read %r has no corresponding label' %
                                 read.id)
            trueLabels.append(trueLabel)
            offsetDeltas = Counter()
            scannedRead = database.scan(read)
            for landmark, trigPoint in database.getScannedPairs(scannedRead):
                delta = scale(trigPoint.offset - landmark.offset,
                              database.distanceBase)
                offsetDeltas[delta] += 1
            allOffsetDeltas.append(offsetDeltas)

        nReads = len(reads)

        if nReads == 0:
            raise ValueError('No sequences were found in %r' % fastaFile)

        # Don't check that len(reads) == len(labels). I.e., extra labels
        # are ignored, to make using this class interactively more convenient.

        # Create an affinity matrix. We initially set all values to 1.0 so
        # we don't need to later initialize the diagonal.
        affinity = np.ones((nReads, nReads))

        for row, offsetDeltas in enumerate(allOffsetDeltas):
            for col in xrange(row + 1, nReads):
                affinity[row, col] = affinity[col, row] = (
                    self.affinityFromOffsetDeltas(
                        allOffsetDeltas[row], allOffsetDeltas[col]))

        self.nReads = nReads
        self.affinity = affinity
        self.trueLabels = trueLabels

    @staticmethod
    def affinityFromOffsetDeltas(a, b):
        """
        Calculate the [0.0, 1.0] affinity between two offset delta Counters.

        @param a: A C{Counter} instance containing counts for the offset deltas
            for a read.
        @param b: A C{Counter} instance containing counts for the offset deltas
            for a read.
        @return: A C{float} from 0.0 to 1.0 (inclusive) of the affinity
            (similarity) of the contents of the two counters.
        """
        # The affinity is twice the size of the intersection divided by the
        # sum of the sizes of the two counters. So if both counters are
        # identical, the affinity will be 1.0 and if they have nothing in
        # common it will be 0.0.
        denom = sum(a.values()) + sum(b.values())
        if denom:
            return 2.0 * sum((a & b).values()) / denom
        else:
            # Neither sequence had any features, so let's say they have no
            # affinity.
            return 0.0


class AffinityPropagationAnalysis(ClusterAnalysis):
    """
    Run affinity propagation cluster analysis.
    """

    def cluster(self):
        """
        Interface to affinity propagation clustering.

        @return: An C{sklearn.cluster.AffinityPropagation} instance.
        """
        affinityPropagation = AffinityPropagation(
            affinity='precomputed').fit(self.affinity)
        self.nClusters = len(affinityPropagation.cluster_centers_indices_)
        self.clusterLabels = affinityPropagation.labels_
        self.affinityPropagation = affinityPropagation
        return affinityPropagation

    def print_(self, fp=sys.stdout):
        """
        Print details of the clustering.

        @param fp: A file pointer to write to.
        """
        try:
            self.clusterLabels
        except AttributeError:
            # Looks like clustering has not been run. Run it.
            self.cluster()

        trueLabels, clusterLabels = self.trueLabels, self.clusterLabels

        print >>fp, 'Estimated number of clusters: %d' % self.nClusters

        print >>fp, 'Homogeneity: %0.3f' % (
            homogeneity_score(trueLabels, clusterLabels))

        print >>fp, 'Completeness: %0.3f' % (
            completeness_score(trueLabels, clusterLabels))

        print >>fp, 'V-measure: %0.3f' % (
            v_measure_score(trueLabels, clusterLabels))

        print >>fp, 'Adjusted Rand Index: %0.3f' % (
            adjusted_rand_score(trueLabels, clusterLabels))

        print >>fp, 'Adjusted Mutual Information: %0.3f' % (
            adjusted_mutual_info_score(trueLabels, clusterLabels))

        print >>fp, 'Silhouette Coefficient: %0.3f' % silhouette_score(
            self.affinity, clusterLabels, metric='sqeuclidean')


class KMeansAnalysis(ClusterAnalysis):
    """
    Interface to k-means clustering.
    """

    def cluster(self, k):
        """
        Cluster the data using the k-means method.

        @return: An C{sklearn.cluster.KMeans} instance.
        """
        self.kmeans = KMeans(n_clusters=k, init='k-means++').fit(self.affinity)
        self.clusterLabels = self.kmeans.labels_
        return self.kmeans

    def print_(self, fp=sys.stdout):
        """
        Print details of the clustering.

        @param fp: A file pointer to write to.
        """
        try:
            trueLabels, clusterLabels = self.trueLabels, self.clusterLabels
        except AttributeError:
            # Looks like clustering has not been run. Ask the user to run it,
            # seeing as we don't want to silently use a default value of k.
            raise RuntimeError('Did you forget to run cluster()?')

        print >>fp, 'Homogeneity: %0.3f' % (
            homogeneity_score(trueLabels, clusterLabels))

        print >>fp, 'Completeness: %0.3f' % (
            completeness_score(trueLabels, clusterLabels))

        print >>fp, 'V-measure: %0.3f' % (
            v_measure_score(trueLabels, clusterLabels))

        print >>fp, 'Adjusted Rand Index: %0.3f' % (
            adjusted_rand_score(trueLabels, clusterLabels))

        print >>fp, 'Adjusted Mutual Information: %0.3f' % (
            adjusted_mutual_info_score(trueLabels, clusterLabels))

        print >>fp, 'Silhouette Coefficient: %0.3f' % silhouette_score(
            self.affinity, clusterLabels, metric='sqeuclidean',
            sample_size=300)
