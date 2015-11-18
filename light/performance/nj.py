from collections import Counter
import numpy as np

from skbio.tree import nj
from skbio import DistanceMatrix, TreeNode

from dark.reads import AARead
from dark.fasta import FastaReads

from light.performance.affinity import affinityMatrix
from light.parameters import FindParameters


class NJTree:

    def __init__(self):
        # You probably don't want to use NJTree() directly. Rather, make an
        # instance by using either NJTree.fromSequences or
        # NJTree.fromDistanceMatrix.
        self.sequences = self.distance = self.tree = self.labels = None

    @classmethod
    def fromSequences(cls, labels, sequences, findParams=None, **kwargs):
        """
        Construct an NJTree instance from some seqeunces.

        @param cls: Our class.
        @param labels: An iterable producing C{str} labels for the sequences.
        @param sequences: Either A C{str} filename of sequences to consider or
            a C{light.reads.Reads} instance.
        @param findParams: An instance of C{FindParameters}.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        if isinstance(sequences, str):
            sequences = FastaReads(sequences, readClass=AARead, upperCase=True)

        new = cls()
        new.sequences = list(sequences)
        new.labels = labels
        findParams = findParams or FindParameters()
        affinity = affinityMatrix(new.sequences, findParams=findParams,
                                  **kwargs)
        new.distance = np.ones(affinity.shape) - affinity
        new.tree = nj(DistanceMatrix(new.distance, labels))
        return new

    @classmethod
    def fromDistanceMatrix(cls, labels, distance, findParams=None, **kwargs):
        """
        Construct an NJTree instance, given a distance matrix.

        @param cls: Our class.
        @param labels: An iterable producing C{str} labels corresponding to the
            rows (equivalently, columns) of the distance matrix.
        @param distance: A square matrix of numeric distances.
        @param findParams: An instance of C{FindParameters}.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        new = cls()
        new.labels = labels
        findParams = findParams or FindParameters()
        new.distance = distance
        new.tree = nj(DistanceMatrix(distance, labels))
        return new

    def canonicalize(self):
        """
        Canonicalize a tree by sorting the children of each internal node so
        they are ordered by 1) node length (the distance to its parent), 2)
        the number of tips they lead to 3) the name of the node (usually only
        tips have names), and 4) the sorted set of names of all tips they lead
        to.

        @return: A new C{NJTree} instance with iternal nodes and children
            ordered canonically.
        """
        def _key(node):
            """
            Compute a canonical sort value for a node.

            @return: A tuple containing the values to sort on (see docstring
                for C{canonicalize}, above).
            """
            return (node.length, node.count(tips=True),
                    node.name or '', sorted(node.subset()))

        def _canonicalize(node):
            """
            Canonicalize node. See docstring for C{canonicalize}, above.
            """
            # This is very inefficient. The key function (above) computes
            # sorted lists of node tip names repeatedly in a naive way.
            # That could be done more efficiently by working from the
            # leaves to the root, combining tip names. For now just use
            # this slow but simple approach.
            if node.is_tip():
                return node.copy()
            else:
                children = list(map(_canonicalize, node.children))
                children.sort(key=_key)
                return TreeNode(children=children, length=node.length,
                                name=node.name)

        new = NJTree()
        new.labels = self.labels
        new.sequences = self.sequences
        new.distance = self.distance
        new.tree = _canonicalize(self.tree)
        return new

    def countClades(self):
        """
        Count all the clades in our tree.

        @return: A C{Counter} instance, whose keys are C{frozenset} instances
            of tip names, corresponding to the clades found under the internal
            nodes of the tree.
        """
        counts = Counter()
        # This is very inefficient (depending on the implementation of the
        # skbio tree), as it recomputes the subset for each node. It would
        # be more efficient to work upwards from each tip, combining
        # subsets as we head to the root.
        for node in self.tree.non_tips(include_self=True):
            counts[node.subset()] += 1
        return counts
