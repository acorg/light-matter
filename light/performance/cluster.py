from collections import Counter
import numpy as np

from sklearn.cluster import AffinityPropagation, KMeans
from sklearn.metrics import (
    homogeneity_score, completeness_score, v_measure_score,
    adjusted_rand_score, adjusted_mutual_info_score, silhouette_score,
    confusion_matrix)

from dark.reads import AARead
from dark.fasta import FastaReads

from light.distance import scale
from light.database import DatabaseSpecifier
from light.backend import Backend
from light.string import MultilineString


class ClusterAnalysis(object):

    def __init__(self, sequences, labels, defaultLabel=None, **kwargs):
        """
        Base class for using cluster analysis to evaluate how well various
        feature finders and database parameter settings can separate a set of
        sequences. The clustering is based on feature offset deltas.

        @param sequences: Either A C{str} filename of sequences to consider or
            a C{light.reads.Reads} instance.
        @param labels: A C{dict} with a label for each sequence id in
            C{sequences}. These are the known categories of each sequence.
        @param defaultLabel: If not C{None}, a label to use for reads whose ids
            are not present in C{labels}. If C{None} and a read id has no label
            a ValueError is raised.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        @raises ValueError: If the id of a read is not in labels and no default
            label has been set, or if there are no reads in C{sequences}.
        """
        if isinstance(sequences, str):
            reads = FastaReads(sequences, readClass=AARead)
        else:
            reads = sequences
        database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(database.params)
        allOffsetDeltas = []
        trueLabels = []

        for read in reads:
            trueLabel = labels.get(read.id, defaultLabel)
            if trueLabel is None:
                raise ValueError('Read %r has no corresponding label' %
                                 read.id)
            trueLabels.append(trueLabel)
            offsetDeltas = Counter()
            scannedRead = backend.scan(read)
            for landmark, trigPoint in backend.getScannedPairs(scannedRead):
                delta = scale(trigPoint.offset - landmark.offset,
                              database.params.distanceBase)
                offsetDeltas[delta] += 1
            allOffsetDeltas.append(offsetDeltas)

        nReads = len(reads)

        if nReads == 0:
            raise ValueError('No sequences were found in %r' % sequences)

        # Don't check that len(reads) == len(labels). I.e., ignore extra labels
        # to make using this class interactively more convenient.

        # Create an affinity matrix. Initially set all values to 1.0 so we
        # don't need to later initialize the diagonal.
        affinity = np.ones((nReads, nReads))

        for row, offsetDeltas in enumerate(allOffsetDeltas):
            for col in range(row + 1, nReads):
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
        # The affinity is the size of the intersection divided by the size
        # of the union of the two counters. If both counters are identical,
        # the affinity will be 1.0 and if they have nothing in common it
        # will be 0.0.
        denom = sum((a | b).values())
        if denom:
            return sum((a & b).values()) / float(denom)
        else:
            # Neither sequence had any features (which means there are no
            # deltas at all). Let's define this as zero affinity.
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
        self.confusionMatrix = confusion_matrix(self.trueLabels,
                                                self.clusterLabels)
        return affinityPropagation

    def print_(self, margin=''):
        """
        Print details of the clustering.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the statistical measures from
            running the affinity propagation clustering.
        """
        try:
            self.clusterLabels
        except AttributeError:
            # Looks like clustering has not been run. Run it.
            self.cluster()

        result = MultilineString(margin=margin)
        append = result.append

        trueLabels, clusterLabels = self.trueLabels, self.clusterLabels

        append('Estimated number of clusters: %d' % self.nClusters)

        append('Homogeneity: %0.3f' % (
            homogeneity_score(trueLabels, clusterLabels)))

        append('Completeness: %0.3f' % (
            completeness_score(trueLabels, clusterLabels)))

        append('V-measure: %0.3f' % (
            v_measure_score(trueLabels, clusterLabels)))

        append('Adjusted Rand Index: %0.3f' % (
            adjusted_rand_score(trueLabels, clusterLabels)))

        append('Adjusted Mutual Information: %0.3f' % (
            adjusted_mutual_info_score(trueLabels, clusterLabels)))

        append('Silhouette Coefficient: %0.3f' % silhouette_score(
            self.affinity, clusterLabels, metric='sqeuclidean'))

        return str(result)


class KMeansAnalysis(ClusterAnalysis):
    """
    Interface to k-means clustering.
    """

    def cluster(self, k):
        """
        Cluster the data using the k-means method.

        @return: An C{sklearn.cluster.KMeans} instance.
        """
        self.kMeans = KMeans(n_clusters=k, init='k-means++').fit(self.affinity)
        self.nClusters = k
        self.clusterLabels = self.kMeans.labels_
        self.confusionMatrix = confusion_matrix(self.trueLabels,
                                                self.clusterLabels)
        return self.kMeans

    def print_(self, margin=''):
        """
        Print details of the clustering.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the statistical measures from
            running the k-means clustering.
        """
        try:
            trueLabels, clusterLabels = self.trueLabels, self.clusterLabels
        except AttributeError:
            # Looks like clustering has not been run. Ask the user to run it,
            # seeing as we don't want to silently use a default value of k.
            raise RuntimeError('Did you forget to run cluster()?')

        result = MultilineString()
        append = result.append

        append('Homogeneity: %0.3f' % (
            homogeneity_score(trueLabels, clusterLabels)))

        append('Completeness: %0.3f' % (
            completeness_score(trueLabels, clusterLabels)))

        append('V-measure: %0.3f' % (
            v_measure_score(trueLabels, clusterLabels)))

        append('Adjusted Rand Index: %0.3f' % (
            adjusted_rand_score(trueLabels, clusterLabels)))

        append('Adjusted Mutual Information: %0.3f' % (
            adjusted_mutual_info_score(trueLabels, clusterLabels)))

        append('Silhouette Coefficient: %0.3f' % silhouette_score(
            self.affinity, clusterLabels, metric='sqeuclidean',
            sample_size=300))

        return str(result)
