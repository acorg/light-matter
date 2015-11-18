from unittest import TestCase
import numpy as np

from light.performance.ultrametricity import ultrametric


class TestUltrametricity(TestCase):
    """
    Tests for the light.performance.ultrametricity function.
    """
    def testYieldsAllPermutations(self):
        """
        If a matrix isn't ultrametric, all triplets must be returned.
        """
        matrix = np.array([[1, 1, 1, 1],
                           [1, 1, 1, 1],
                           [1, 1, 1, 1],
                           [1, 1, 1, 1]])
        result = len(list(ultrametric(matrix)))
        self.assertEqual(24, result)

    def testYieldsAllNonUltrametric(self):
        """
        All non-ultrametric triplets must be returned.
        """
        matrix = np.array([[1, 2, 1, 1],
                           [1, 1, 1, 1],
                           [1, 1, 1, 1],
                           [1, 1, 1, 1]])
        result = len(list(ultrametric(matrix)))
        self.assertEqual(22, result)
