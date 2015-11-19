from unittest import TestCase
import numpy as np

from light.performance.ultrametricity import ultrametric


class TestUltrametricity(TestCase):
    """
    Tests for the light.performance.ultrametricity function.
    """
    def testYieldsNoPermutations(self):
        """
        If a matrix is ultrametric, no triplets must be returned.
        """
        matrix = np.array([[0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0]])
        result = list(ultrametric(matrix))
        self.assertEqual([], result)

    def testYieldsAllNonUltrametric(self):
        """
        All non-ultrametric triplets must be returned.
        """
        matrix = np.array([[0, 1, 2],
                           [1, 0, 1],
                           [2, 1, 0]])
        result = list(ultrametric(matrix))
        self.assertEqual([(0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1)], result)
