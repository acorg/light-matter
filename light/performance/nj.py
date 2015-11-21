from collections import Counter
from io import StringIO

import matplotlib.pyplot as plt
import numpy as np

from Bio import Phylo

from skbio.tree import nj, majority_rule
from skbio import DistanceMatrix, TreeNode

from dark.reads import AARead
from dark.fasta import FastaReads

from light.performance.affinity import affinityMatrix
from light.parameters import FindParameters


# The default standard deviataion for use in perturbDistanceMatrix.
_DEFAULT_STDDEV = 0.05


def perturbDistanceMatrix(matrix, stddev=None):
    """
    Add a small amount of symmetric off-diagonal noise to a square distance
    matrix.

    @param matrix: either a square matrix of numeric values or a square
        C{np.array}.
    @param stddev: The C{float} standard deviation of the noise to add to
        off-diagonal elements. If C{None}, use _DEFAULT_STDDEV.
    @return: An C{np.array} containing the same diagonal values as C{matrix}
        and with some normally distribued noise added (symmetrically) to the
        off-diagonal elements. Off diagonal values in the return result will
        not be less than zero or greater than one.
    """
    stddev = _DEFAULT_STDDEV if stddev is None else stddev
    normal = np.random.normal
    newMatrix = np.copy(matrix)
    n = newMatrix.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            value = newMatrix[i, j] + normal(loc=0.0, scale=stddev)
            if value < 0.0:
                value = 0.0
            elif value > 1.0:
                value = 1.0
            newMatrix[i, j] = newMatrix[j, i] = value
    return newMatrix


class NJTree:

    def __init__(self):
        # You probably don't want to use NJTree() directly. Rather, make an
        # instance by using either NJTree.fromSequences or
        # NJTree.fromDistanceMatrix.
        self.sequences = self.distance = self.tree = self.labels = None
        self.supportIterations = 0

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
        @return: An C{NJTree} instance.
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
    def fromDistanceMatrix(cls, labels, distance):
        """
        Construct an NJTree instance, given a distance matrix.

        @param cls: Our class.
        @param labels: An iterable producing C{str} labels corresponding to the
            rows (equivalently, columns) of the distance matrix.
        @param distance: A square matrix of numeric distances.
        @return: An C{NJTree} instance.
        """
        new = cls()
        new.labels = labels
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
                new = TreeNode(children=children, length=node.length,
                               name=node.name)
                try:
                    new.support = node.support
                except AttributeError:
                    pass
                return new

        new = NJTree()
        new.labels = self.labels
        new.sequences = self.sequences
        new.distance = self.distance
        new.tree = _canonicalize(self.tree)
        new.supportIterations = self.supportIterations
        return new

    def root(self, nodes):
        """
        Root a tree by the lowest common ancestor of a list of nodes.

        @param nodes: Either a C{list} of C{str} tree tip names or a C{list} of
            C{TreeNode} instances.
        @raise MissingNodeError: If a non-existent node name is present in
            C{nodes}.
        @return: A new C{NJTree} instance rooted at lowest common ancestor of
            the specified node.
        """
        lca = self.tree.lca(nodes)
        if lca.is_tip():
            # Can't root at a tip, let's use its parent.
            lca = lca.parent

        if lca.parent is None:
            # The lca is the root of the tree. Nothing to do.
            return self

        new = NJTree()
        new.labels = self.labels
        new.sequences = self.sequences
        new.distance = self.distance
        new.supportIterations = self.supportIterations
        new.tree = self.tree.root_at(lca)
        return new

    def countClades(self):
        """
        Count all the clades in our tree.

        @return: A C{Counter} instance, whose keys are C{frozenset} instances
            of tip names, corresponding to the clades decendant from internal
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

    @staticmethod
    def addSupportToNode(node, support):
        """
        Increment the support level for the clade beneath a node.

        @param node: A C{TreeNode} instance.
        @param support: The numeric support value to add to the existing
            support.
        """
        try:
            node.support += support
        except AttributeError:
            node.support = support

    def supportForNode(self, node):
        """
        Get the support level for a node.

        @param node: A C{TreeNode} instance.
        @return: The C{float} support, ranging from 0.0 to 1.0.
        """
        # Explicitly check if supportIterations is non-zero rather than
        # dividing by it and also catching ZeroDivisionError because
        # node.support is sometimes a numpy.float64 on return from
        # majority_rule (called by self.consensusTrees, below). I don't
        # know why that happens. Anyway, if you divide a numpy.float64 by
        # zero you don't get a ZeroDivisionError, you get 'inf' and a
        # warning.
        if self.supportIterations:
            try:
                return node.support / self.supportIterations
            except AttributeError:
                pass
        return 0.0

    def newick(self):
        """
        Get our tree in Newick format, with clade support in comment fields.

        @return: A Newick format C{str} with clade support in square bracketed
            comment fields if support has been added.
        """

        def _newick(node):
            """
            See docstring for C{newick} above.

            @param node: A C{TreeNode} instance.
            """
            if node.is_tip():
                return '%s:%f' % (node.name, node.length)
            else:
                children = ','.join(_newick(child) for child in node.children)
                support = ('[%d]' % int(self.supportForNode(node) * 100)
                           if self.supportIterations else '')
                suffix = '' if node.length is None else ':%s' % node.length
                return '(%s)%s%s' % (children, support, suffix)

        result = _newick(self.tree)
        if self.supportIterations:
            # Remove the final [100] support which is always present because
            # the root of all topologies has all tips as descendants.
            assert result.endswith('[100]')
            result = result[:-5]

        return result + ';\n'

    def addSupport(self, iterations, stddev=None):
        """
        Make new perturbed distance matrices, create NJ trees from them, and
        add their clade counts to the internal nodes of this tree.

        @param iterations: The C{int} number of times to create new NJ trees
            via a perturbed distance matrix.
        @param stddev: The C{float} standard deviation of the noise to add to
            off-diagonal distance matrix elements before creating a new NJ
            tree to add to the support.
        """
        cladeCounts = Counter()
        for _ in range(iterations):
            distance = perturbDistanceMatrix(self.distance, stddev)
            new = NJTree.fromDistanceMatrix(self.labels, distance)
            cladeCounts.update(new.countClades())

        # Visit each non-tip node and increment its support based on
        # the clade counts just computed.
        for node in self.tree.non_tips(include_self=True):
            self.addSupportToNode(node, cladeCounts[node.subset()])

        self.supportIterations += iterations

    def plot(self, filename=None, show=True, **kwargs):
        """
        Use the Phylo package to plot the NJ tree, showing branch support.

        @param filename: A C{str} file name to which to write the tree image.
        @param show: If C{True} the image will be displayed. This is only
            useful when C{filename} is not C{None}.
        @param kwargs: Additional (optional) arguments to be passed to savefig.
            For available options, see:
            http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.savefig
        """
        fp = StringIO(self.newick())
        Phylo.draw(Phylo.read(fp, 'newick', comments_are_confidence=True),
                   do_show=show)
        if filename:
            plt.gcf().savefig(filename, **kwargs)
        if not show:
            plt.close()

    def consensusTrees(self, iterations, stddev=None):
        """
        Make new perturbed distance matrices, create NJ trees from them, and
        calculate consensus trees based on our original tree and the new ones.

        @param iterations: The C{int} number of times to create a new NJ tree
            via a perturbed distance matrix.
        @param stddev: The C{float} standard deviation of the noise to add to
            off-diagonal distance matrix elements before creating a new NJ
            tree to add to the support.
        @return: A C{list} of new C{NJTree} instances, holding the consensus
            trees produced by C{skbio.tree.majority_rule}.
        """
        trees = [self.tree]
        for _ in range(iterations):
            distance = perturbDistanceMatrix(self.distance, stddev)
            new = NJTree.fromDistanceMatrix(self.labels, distance)
            trees.append(new.tree)

        result = []
        for tree in majority_rule(trees):
            new = NJTree()
            new.labels = self.labels
            new.sequences = self.sequences
            new.distance = self.distance
            new.supportIterations += iterations + 1
            new.tree = tree
            result.append(new)

        return result
