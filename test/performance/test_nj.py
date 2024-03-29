import six
from collections import Counter
from unittest import TestCase
import numpy as np
from six import StringIO

from skbio.tree._exception import MissingNodeError
from skbio.stats.distance import DissimilarityMatrixError
from skbio import TreeNode

from dark.reads import Reads, AARead

from light.performance.nj import NJTree, perturbDistanceMatrix


class TestNJTree(TestCase):
    """
    Tests for the light.performance.nj.NJTree class.
    """

    def testNewTreeHasExpectedAttributes(self):
        """
        A new NJTree instance must have the expected attributes.
        """
        njtree = NJTree()
        self.assertEqual(0, njtree.supportIterations)
        self.assertIs(None, njtree.sequences)
        self.assertIs(None, njtree.distance)
        self.assertIs(None, njtree.tree)
        self.assertIs(None, njtree.labels)

    def testFromSequencesWithNoSequencesRaisesDissimilarityMatrixError(self):
        """
        If an attempt is made to build an NJTree instance from a set of
        sequences that is empty, a DissimilarityMatrixError must be raised.
        """
        sequences = Reads()
        error = '^Data must be at least 1x1 in size\.$'
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromSequences, [], sequences,
                              landmarks=['AlphaHelix'])

    def testFromSequencesWithNoLabels(self):
        """
        If an attempt is made to build an NJTree instance and the number
        of labels does not match the number of sequences, a
        DissimilarityMatrixError must be raised.
        """
        sequences = Reads()
        sequences.add(AARead('id', 'A'))
        error = ('^The number of IDs \(0\) must match the number of '
                 'rows/columns in the data \(1\)\.$')
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromSequences, [], sequences,
                              landmarks=['AlphaHelix'])

    def testFromSequencesWithOneSequenceRaisesValueError(self):
        """
        If an attempt is made to build an NJTree instance from just one
        sequence, a ValueError must be raised.
        """
        sequences = Reads()
        sequences.add(AARead('id', 'A'))
        error = ('^Distance matrix must be at least 3x3 to generate a '
                 'neighbor joining tree\.$')
        six.assertRaisesRegex(self, ValueError, error,
                              NJTree.fromSequences, ['x'], sequences,
                              landmarks=['AlphaHelix'])

    def testFromSequencesWithTwoSequencesRaisesValueError(self):
        """
        If an attempt is made to build an NJTree instance from just two
        sequences, a ValueError must be raised.
        """
        sequences = Reads()
        sequences.add(AARead('id1', 'A'))
        sequences.add(AARead('id2', 'B'))
        error = ('^Distance matrix must be at least 3x3 to generate a '
                 'neighbor joining tree\.$')
        six.assertRaisesRegex(self, ValueError, error,
                              NJTree.fromSequences, ['x', 'y'], sequences,
                              landmarks=['AlphaHelix'])

    def testFromThreeSequences(self):
        """
        If three sequences with no features are used to create an NJTree
        instance, the instance must 1) have a distance matrix that is zero
        on the diagonal and ones elsewhere, 2) save the labels, and 3) produce
        a simple tree with three children.
        """
        sequences = Reads()
        sequences.add(AARead('id1', 'A'))
        sequences.add(AARead('id2', 'A'))
        sequences.add(AARead('id3', 'A'))
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromSequences(labels, sequences,
                                      landmarks=['AlphaHelix'])
        self.assertTrue(np.array_equal(
            [
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 0],
            ],
            njtree.distance))

        self.assertIs(labels, njtree.labels)
        self.assertEqual(['x:0.5;\n', 'y:0.5;\n', 'z:0.5;\n'],
                         sorted(str(child) for child in njtree.tree.children))

    def testFromEmptyDistanceMatrix(self):
        """
        If an attempt is made to build an NJTree instance from an empty
        distance matrix, a DissimilarityMatrixError must be raised.
        """
        error = '^Data must be at least 1x1 in size\.$'
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromDistanceMatrix, [], [])

    def testFromDistanceMatrixOneByOne(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix of size 1x1, an assertion error must be raised.
        """
        error = '^Data must have exactly two dimensions\.$'
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromDistanceMatrix, [], [1.0])

    def testFromDistanceMatrixTwoByTwo(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix that does not have zeroes on the diagonal, a
        ValueError must be raised.
        """
        error = ('^Distance matrix must be at least 3x3 to generate a '
                 'neighbor joining tree\.$')
        six.assertRaisesRegex(self, ValueError, error,
                              NJTree.fromDistanceMatrix,
                              ['a', 'b'], [[0, 1], [1, 0]])

    def testFromDistanceMatrixWithNoLabels(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix but no labels are given, a DissimilarityMatrixError must be
        raised.
        """
        error = ('^The number of IDs \(0\) must match the number of '
                 'rows/columns in the data \(3\)\.$')
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromDistanceMatrix,
                              [], [[0, 1, 1], [1, 0, 1], [1, 1, 0]])

    def testFromNonHollowDistanceMatrix(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix that does not have zeroes on the diagonal, a
        DissimilarityMatrixError must be raised.
        """
        error = ('^Data must be hollow \(i.e., the diagonal can only contain '
                 'zeros\)\.$')
        six.assertRaisesRegex(self, DissimilarityMatrixError, error,
                              NJTree.fromDistanceMatrix,
                              [], [[1, 1], [1, 1]])

    def testFromDistanceMatrixThreeByThree(self):
        """
        If an NJTree instance is built from a 3x3 distance matrix, the
        instance must 1) save the distance matrix, 2) save the labels,
        and 3) produce a simple tree with three children, and 4) have no
        sequences.
        """
        distance = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)

        self.assertIs(distance, njtree.distance)
        self.assertIs(labels, njtree.labels)
        self.assertEqual(['x:0.5;\n', 'y:0.5;\n', 'z:0.5;\n'],
                         sorted(str(child) for child in njtree.tree.children))
        self.assertIs(None, njtree.sequences)

    def testWithNoSupportAllNodesHaveSupportOfZero(self):
        """
        If three sequences with no features are used to create an NJTree
        the child nodes in the resulting tree must all have support of 0.0.
        """
        distance = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        self.assertEqual(
            [0.0, 0.0, 0.0],
            [njtree.supportForNode(child) for child in njtree.tree.children])

    def testCanonicalizeByNodeLength(self):
        """
        In forming a canonical tree, child nodes must be sorted by length.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(length=13),
            TreeNode(length=11),
            TreeNode(length=18),
            TreeNode(length=14)])

        self.assertEqual(
            [13, 11, 18, 14],
            [child.length for child in njtree.tree.children])
        self.assertEqual(
            [11, 13, 14, 18],
            [child.length for child in njtree.canonicalize().tree.children])

    def testCanonicalizeByNumberOfTips(self):
        """
        In forming a canonical tree, child nodes must be sorted by number
        of tips (assuming child lengths are all equal).
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(children=[
                TreeNode(),
                TreeNode(),
                TreeNode(),
            ]),
            TreeNode(children=[
                TreeNode(),
                TreeNode(),
                TreeNode(),
                TreeNode(),
                TreeNode(),
            ]),
            TreeNode(children=[
                TreeNode(),
                TreeNode(),
            ]),
        ])

        self.assertEqual(
            [3, 5, 2],
            [len(child.children) for child in njtree.tree.children])
        self.assertEqual(
            [2, 3, 5],
            [len(child.children)
             for child in njtree.canonicalize().tree.children])

    def testCanonicalizeByNodeName(self):
        """
        In forming a canonical tree, child nodes must be sorted by name if node
        lengths and number of tips are equal.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        self.assertEqual(
            ['c', 'd', 'b', 'a'],
            [child.name for child in njtree.tree.children])
        self.assertEqual(
            ['a', 'b', 'c', 'd'],
            [child.name for child in njtree.canonicalize().tree.children])

    def testCanonicalizeByTipSubset(self):
        """
        In forming a canonical tree, child nodes must be sorted by the names of
        the set of tips they lead to, if all else is equal.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(children=[
                TreeNode(name='d'),
                TreeNode(name='e'),
                TreeNode(name='f'),
            ]),
            TreeNode(children=[
                TreeNode(name='g'),
                TreeNode(name='h'),
                TreeNode(name='i'),
            ]),
            TreeNode(children=[
                TreeNode(name='a'),
                TreeNode(name='b'),
                TreeNode(name='c'),
            ]),
        ])

        self.assertEqual(
            ['d', 'e', 'f', 'g', 'h', 'i', 'a', 'b', 'c'],
            [grandchild.name for child in njtree.tree.children
             for grandchild in child.children])
        self.assertEqual(
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'],
            [grandchild.name for child in njtree.canonicalize().tree.children
             for grandchild in child.children])

    def testCanonicalizeByLengthAheadOfNumberOfTips(self):
        """
        In forming a canonical tree, child nodes must be sorted by length
        in preference to number of tips.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(
                length=13,
                children=[
                    TreeNode(),
                    TreeNode(),
                    TreeNode(),
                ],
            ),
            TreeNode(
                length=11,
                children=[
                    TreeNode(),
                    TreeNode(),
                    TreeNode(),
                    TreeNode(),
                ],
            ),
            TreeNode(
                length=12,
                children=[
                    TreeNode(),
                    TreeNode(),
                ],
            ),
        ])

        self.assertEqual(
            [13, 11, 12],
            [child.length for child in njtree.tree.children])
        self.assertEqual(
            [11, 12, 13],
            [child.length for child in njtree.canonicalize().tree.children])

    def testCanonicalizeByNumberOfTipsAheadOfName(self):
        """
        In forming a canonical tree, child nodes must be sorted by number of
        tips in preference to name.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='a',
                     children=[
                         TreeNode(),
                         TreeNode(),
                         TreeNode(),
                     ]),
            TreeNode(name='b',
                     children=[
                         TreeNode(),
                     ]),
            TreeNode(name='c',
                     children=[
                         TreeNode(),
                         TreeNode(),
                     ]),
        ])

        self.assertEqual(
            ['a', 'b', 'c'],
            [child.name for child in njtree.tree.children])
        self.assertEqual(
            ['b', 'c', 'a'],
            [child.name for child in njtree.canonicalize().tree.children])

    def testCanonicalizeByNameAheadOfNamesOfDescendants(self):
        """
        In forming a canonical tree, child nodes must be sorted by name in
        preference to the sorted names of all their descendants.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='b',
                     children=[
                         TreeNode(name='e'),
                         TreeNode(name='f'),
                         TreeNode(name='g'),
                     ]),
            TreeNode(name='c',
                     children=[
                         TreeNode(name='h'),
                         TreeNode(name='i'),
                         TreeNode(name='j'),
                     ]),
            TreeNode(name='a',
                     children=[
                         TreeNode(name='k'),
                         TreeNode(name='l'),
                         TreeNode(name='m'),
                     ]),
        ])

        self.assertEqual(
            ['b', 'c', 'a'],
            [child.name for child in njtree.tree.children])
        self.assertEqual(
            ['a', 'b', 'c'],
            [child.name for child in njtree.canonicalize().tree.children])

    def testRootByMultipleNodeNames(self):
        """
        Rooting by multiple node names must work.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        self.assertEqual(
            ['c', 'd', 'b', 'a'],
            [child.name for child in njtree.root(['a', 'b']).tree.children])

    def testRootByOneNodeName(self):
        """
        Rooting by one node name must work.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        self.assertEqual(
            ['c', 'd', 'b', 'a'],
            [child.name for child in njtree.root(['a']).tree.children])

    def testRootByInexistentNodeNameMustRaiseError(self):
        """
        Rooting by an inexistent node name must raise an exception.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        error = 'Node f is not in self'

        six.assertRaisesRegex(self, MissingNodeError, error, njtree.root,
                              ['f'])

    def testRootByOneTreeNode(self):
        """
        Rooting by one TreeNode must work.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        node = njtree.tree.find('a')

        self.assertEqual(
            ['c', 'd', 'b', 'a'],
            [child.name for child in njtree.root([node]).tree.children])

    def testRootByTwoTreeNodes(self):
        """
        Rooting by two TreeNodes must work.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='c'),
            TreeNode(name='d'),
            TreeNode(name='b'),
            TreeNode(name='a')])

        node1 = njtree.tree.find('a')
        node2 = njtree.tree.find('b')

        self.assertEqual(
            ['c', 'd', 'b', 'a'],
            [child.name for child in
             njtree.root([node1, node2]).tree.children])

    def testCountCladesEmptyTree(self):
        """
        In a tree with no children, there are no clades.
        """
        njtree = NJTree()
        njtree.tree = TreeNode()
        self.assertEqual(Counter(), njtree.countClades())

    def testCountCladesOneChild(self):
        """
        In a tree with one child, there is one clade.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(name='a'),
        ])
        self.assertEqual(
            {
                frozenset(['a']): 1,
            },
            njtree.countClades()
        )

    def testCountCladesTwoChildren(self):
        """
        In a tree with two children, one of which has two children, there are
        two clades.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(children=[
                TreeNode(name='a'),
                TreeNode(name='b'),
            ]),
            TreeNode(name='c'),
        ])
        self.assertEqual(
            {
                frozenset(['a', 'b']): 1,
                frozenset(['a', 'b', 'c']): 1,
            },
            njtree.countClades()
        )

    def testCountCladesTwoChildrenOneGrandchild(self):
        """
        In a tree with two children, one of which has two children (one with a
        grandchild), there are three clades.
        """
        njtree = NJTree()
        njtree.tree = TreeNode(children=[
            TreeNode(children=[
                TreeNode(name='a'),
                TreeNode(children=[
                    TreeNode(name='b'),
                    TreeNode(name='c'),
                ]),
            ]),
            TreeNode(name='d'),
        ])
        self.assertEqual(
            {
                frozenset(['b', 'c']): 1,
                frozenset(['a', 'b', 'c']): 1,
                frozenset(['a', 'b', 'c', 'd']): 1,
            },
            njtree.countClades()
        )

    def testDictCanonicalizes(self):
        """
        When hashing NJTrees, the hash function must canonicalize the trees.
        So, adding the same tree (modulo canonicalization) twice to a dict must
        result in a dict of size one.
        """
        njtree1 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((b:0.4,a:0.9):0.7);'))
        njtree2 = NJTree()
        njtree2.tree = TreeNode.read(StringIO('((a:0.9,b:0.4):0.7);'))
        self.assertEqual(1, len(dict.fromkeys([njtree1, njtree2])))

    def testSetCanonicalizes(self):
        """
        When putting NJTrees into sets, the hash function must canonicalize
        the trees. So, adding the same tree (modulo canonicalization) twice to
        a set must result in a set of size one.
        """
        njtree1 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a:1,b:2)c);'))
        njtree2 = NJTree()
        njtree2.tree = TreeNode.read(StringIO('((b:2,a:1)c);'))
        self.assertEqual(1, len({njtree1, njtree2}))

    def testSupportForNodeWhenNoSuppportAdded(self):
        """
        If no support has been added to a tree, supportForNode must return
        zero for all nodes.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        self.assertEqual(0, njtree.supportForNode(njtree.tree))
        self.assertEqual([0, 0, 0],
                         [njtree.supportForNode(child)
                          for child in njtree.tree.children])

    def testAddSupportIncrementsSupportIterations(self):
        """
        When support has been added to a tree, its supportIterations attribute
        must be incremented correctly.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        njtree.addSupport(2)
        self.assertEqual(2, njtree.supportIterations)

    def testSupportForNodeIsOneAtRoot(self):
        """
        When support has been added to a tree, the root node must have support
        of 1.0. This is because all trees will have all tips under their root,
        regardless of their topologies.
        """
        distance = [
            [0.0, 0.5, 0.4],
            [0.5, 0.0, 0.1],
            [0.4, 0.1, 0.0],
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        njtree.addSupport(10)
        self.assertEqual(1.0, njtree.supportForNode(njtree.tree))

    def testNewickWithNoSuppportAdded(self):
        """
        If no support has been added to a tree, newick must return the expected
        string (with no support values).
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        # The order in the Newick string seems deterministic, and according
        # to the skbio docs this is the case.
        self.assertEqual('(y:0.500000,x:0.500000,z:0.500000);\n',
                         njtree.newick())

    def testNewickWithSupportAdded(self):
        """
        If support has been added to a tree, the Newick string must contain
        square brackets (Newick comment fields in which support is shown).
        """
        distance = [
            [0, 5, 5, 5],
            [5, 0, 5, 1],
            [5, 5, 0, 5],
            [5, 1, 5, 0],
        ]
        labels = ['w', 'x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        njtree.addSupport(2)
        newick = njtree.newick(includeSupport=True)
        self.assertNotEqual(-1, newick.find('['))
        self.assertNotEqual(-1, newick.find(']'))

    def testNewickWithNoSupportIncludedWhenSupportHasBeenAdded(self):
        """
        If support has been added to a tree, but includeSupport is passed to
        newick() as False, the Newick string must not contain any square
        brackets (Newick comment fields in which support is shown).
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        njtree.addSupport(2)
        newick = njtree.newick(includeSupport=False)
        self.assertEqual(-1, newick.find('['))
        self.assertEqual(-1, newick.find(']'))

    def testOnlyOneConsensusTreeWithZeroIterations(self):
        """
        When consensusTrees is passed an iteration count of zero, only one
        consensus tree must be returned.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        consensusTrees = njtree.consensusTrees(0)
        self.assertEqual(1, len(consensusTrees))

    def testConsensusTreeWithZeroIterationsHasSupportOneForAllChildren(self):
        """
        When consensusTrees is passed an iteration count of zero, all children
        in the consensus tree must have support of zero.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        consensusTrees = njtree.consensusTrees(0)
        self.assertEqual([0, 0, 0],
                         [njtree.supportForNode(child)
                          for child in consensusTrees[0].tree.children])

    def testRobinsonFoulds(self):
        """
        The correct Robinson Foulds distance must be calculated.
        """
        njtree1 = NJTree()
        njtree2 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        njtree2.tree = TreeNode.read(StringIO('(((a,b),c),d);'))
        distance = njtree1.robinsonFoulds(njtree2)
        self.assertEqual(2.0, distance)

    def testRobinsonFouldsReversed(self):
        """
        The correct Robinson Foulds distance must be calculated if the trees
        are passed in reverse order.
        """
        njtree1 = NJTree()
        njtree2 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        njtree2.tree = TreeNode.read(StringIO('(((a,b),c),d);'))
        distance = njtree2.robinsonFoulds(njtree1)
        self.assertEqual(2.0, distance)

    def testRobinsonFouldsProportionTrue(self):
        """
        The correct Robinson Foulds distance must be calculated with
        proportion=True.
        """
        njtree1 = NJTree()
        njtree2 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        njtree2.tree = TreeNode.read(StringIO('(((a,b),c),d);'))
        distance = njtree1.robinsonFoulds(njtree2, proportion=True)
        self.assertEqual(0.5, distance)

    def testRobinsonFouldsProportionTrueReversed(self):
        """
        The correct Robinson Foulds distance must be calculated with
        proportion=True, if the trees are passed in reverse order.
        """
        njtree1 = NJTree()
        njtree2 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        njtree2.tree = TreeNode.read(StringIO('(((a,b),c),d);'))
        distance = njtree2.robinsonFoulds(njtree1, proportion=True)
        self.assertEqual(0.5, distance)

    def testRobinsonFouldsCompareAgainstItself(self):
        """
        If a tree is compared against itself, the Robinson Foulds distance must
        be 0.0.
        """
        njtree1 = NJTree()
        njtree2 = NJTree()
        njtree1.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        njtree2.tree = TreeNode.read(StringIO('((a,b),(c,d));'))
        distance = njtree1.robinsonFoulds(njtree2)
        self.assertEqual(0.0, distance)

    def testGenerateTreesZero(self):
        """
        When generateTrees is called with a zero argument, it should return an
        empty list.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        self.assertEqual([], list(njtree.generateTrees(0)))

    def testGenerateTreesMultiple(self):
        """
        When generateTrees is called with a non-zero argument, it should
        return a list with the expected number of items, which must all
        be instances of NJTree.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        result = list(njtree.generateTrees(5))
        self.assertEqual(5, len(result))
        self. assertTrue(all(isinstance(t, NJTree) for t in result))

    def testTreeSampleZero(self):
        """
        When treeSample is called with a zero argument, it should return an
        empty list.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        self.assertEqual([], njtree.treeSample(0))

    def testTreeSampleReturnType(self):
        """
        The treeSample method must return a list of tuples, each containing
        an NJTree and a count.
        """
        distance = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0]
        ]
        labels = ['x', 'y', 'z']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        result = njtree.treeSample(10)
        for tree, count in result:
            self.assertTrue(isinstance(tree, NJTree))
            self.assertTrue(isinstance(count, int))

    def testTreeSampleSortedCounts(self):
        """
        The treeSample method returns a list of tuples, the second element of
        which is a count. The counts must be sorted from highest to lowest.
        """
        distance = [
            [0, 1, 1, 1, 1],
            [1, 0, 1, 1, 1],
            [1, 1, 0, 1, 1],
            [1, 1, 1, 0, 1],
            [1, 1, 1, 1, 0],
        ]
        labels = ['a', 'b', 'c', 'd', 'e']
        njtree = NJTree.fromDistanceMatrix(labels, distance)
        # This test is not 100% reliable because there is a tiny chance the
        # treeSample method will repeatedly generate the same tree, even
        # though we make 50 trees and pass a high (0.3) standard deviation.
        # It is also possible that the counts come back sorted just by
        # chance even if the treeSample method does not sort them (in which
        # case the test will pass even though the treeSample code is wrong).
        #
        # To make these errors even more unlikely, I have put in an assert
        # on the number of trees that come back and I do the test multiple
        # times. If the assert fails it doesn't mean that treeSample is
        # broken, but it's highly likely that it is! It's also extremely
        # unlikely that the counts could be sorted correctly by chance on
        # every run. So if the test passes it doesn't mean that treeSample
        # isn't broken, but it's highly unlikely that it is! Phew.
        for _ in range(10):
            counts = [count for (tree, count) in njtree.treeSample(50, 0.3)]
            self.assertTrue(len(counts) > 2)
            self.assertEqual(sorted(counts, reverse=True), counts)


class TestPerturbDistanceMatrix(TestCase):
    """
    Tests for the light.performance.nj.perturbDistanceMatrix function.
    """

    def testZeroStandardDeviation(self):
        """
        If zero is passed to perturbDistanceMatrix as the standard deviation
        of the noise, ValueError must be raised.
        """
        distance = [
            [0.0, 0.7, 0.3],
            [0.7, 0.0, 0.1],
            [0.3, 0.1, 0.0],
        ]
        error = '^scale <= 0$'
        six.assertRaisesRegex(self, ValueError, error, perturbDistanceMatrix,
                              distance, 0.0)

    def testZeroStandardDeviationx(self):
        """
        If zero is passed to perturbDistanceMatrix as the standard deviation
        of the noise, ValueError must be raised.
        """
        distance = [
            [0.0, 0.7],
            [0.7, 0.0],
        ]
        error = '^scale <= 0$'
        six.assertRaisesRegex(self, ValueError, error, perturbDistanceMatrix,
                              distance, 0.0)

    def testResultTypeAndShape(self):
        """
        perturbDistanceMatrix must return an np.ndarray of same shape as
        its matrix argument.
        """
        distance = np.array([
            [0.0, 0.7, 0.1],
            [0.7, 0.0, 0.4],
            [0.1, 0.4, 0.0],
        ])
        result = perturbDistanceMatrix(distance)
        self.assertEqual((3, 3), result.shape)
        self.assertTrue(isinstance(result, np.ndarray))

    def testDiagonalUntouched(self):
        """
        perturbDistanceMatrix must not change any values on the diagonal.
        """
        # Use non-zero diagonal values to get more confidence that those
        # values are not being altered.
        distance = np.array([
            [0.2, 0.7, 0.1],
            [0.7, 0.5, 0.4],
            [0.1, 0.4, 0.9],
        ])
        result = perturbDistanceMatrix(distance, 0.3)
        self.assertEqual([0.2, 0.5, 0.9], [result[i][i] for i in range(3)])
