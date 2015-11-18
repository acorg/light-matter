from collections import Counter
from unittest import TestCase
import numpy as np

from skbio.stats.distance import DissimilarityMatrixError
from skbio import TreeNode

from dark.reads import Reads, AARead

from light.performance.nj import NJTree


class TestNJTree(TestCase):
    """
    Tests for the light.performance.nj.NJTree class.
    """

    def testFromSequencesWithNoSequencesRaisesDissimilarityMatrixError(self):
        """
        If an attempt is made to build an NJTree instance from a set of
        sequences that is empty, a DissimilarityMatrixError must be raised.
        """
        sequences = Reads()
        error = '^Data must be at least 1x1 in size\.$'
        self.assertRaisesRegex(DissimilarityMatrixError, error,
                               NJTree.fromSequences, [], sequences,
                               landmarkNames=['AlphaHelix'])

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
        self.assertRaisesRegex(DissimilarityMatrixError, error,
                               NJTree.fromSequences, [], sequences,
                               landmarkNames=['AlphaHelix'])

    def testFromSequencesWithOneSequenceRaisesValueError(self):
        """
        If an attempt is made to build an NJTree instance from just one
        sequence, a ValueError must be raised.
        """
        sequences = Reads()
        sequences.add(AARead('id', 'A'))
        error = ('^Distance matrix must be at least 3x3 to generate a '
                 'neighbor joining tree\.$')
        self.assertRaisesRegex(ValueError, error,
                               NJTree.fromSequences, ['x'], sequences,
                               landmarkNames=['AlphaHelix'])

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
        self.assertRaisesRegex(ValueError, error,
                               NJTree.fromSequences, ['x', 'y'], sequences,
                               landmarkNames=['AlphaHelix'])

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
                                      landmarkNames=['AlphaHelix'])
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
        self.assertRaisesRegex(DissimilarityMatrixError, error,
                               NJTree.fromDistanceMatrix, [], [])

    def testFromDistanceMatrixOneByOne(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix of size 1x1, an assertion error must be raised.
        """
        error = '^Data must have exactly two dimensions\.$'
        self.assertRaisesRegex(DissimilarityMatrixError, error,
                               NJTree.fromDistanceMatrix, [], [1.0])

    def testFromDistanceMatrixTwoByTwo(self):
        """
        If an attempt is made to build an NJTree instance from a distance
        matrix that does not have zeroes on the diagonal, a
        ValueError must be raised.
        """
        error = ('^Distance matrix must be at least 3x3 to generate a '
                 'neighbor joining tree\.$')
        self.assertRaisesRegex(ValueError, error,
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
        self.assertRaisesRegex(DissimilarityMatrixError, error,
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
        self.assertRaisesRegex(DissimilarityMatrixError, error,
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
