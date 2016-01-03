from unittest import TestCase

from light.performance.overlap import SSAARead

from light.features import Landmark
from light.landmarks.predicted_structure import PredictedStructure


class TestGetStructureOffsets(TestCase):
    """
    Tests for the Landmark.PredictedStructure.getStructureOffsets function.
    """
    def testNoStructureSequence(self):
        """
        An empty structure sequence must return an empty dict.
        """
        landmark = PredictedStructure()
        result = landmark.getStructureOffsets('')
        self.assertEqual({}, result)

    def testIdenticalStructure(self):
        """
        A structure sequence that only consists of one type of structure must
        lead to the correct dict.
        """
        landmark = PredictedStructure()
        result = landmark.getStructureOffsets('GGGGGGGGG')
        self.assertEqual({'G': [[0, 8]]}, result)

    def testFirstItemDifferent(self):
        """
        If the first item in the structure sequence is different from the rest,
        the right dict must be returned.
        """
        landmark = PredictedStructure()
        result = landmark.getStructureOffsets('HGGGGGGGG')
        self.assertEqual({'G': [[1, 8]],
                          'H': [[0, 0]]}, result)

    def testLastItemDifferent(self):
        """
        If the last item in the structure sequence is different from the rest,
        the right dict must be returned.
        """
        landmark = PredictedStructure()
        result = landmark.getStructureOffsets('GGGGGGGGH')
        self.assertEqual({'G': [[0, 7]],
                          'H': [[8, 8]]}, result)

    def testLongStructureSequence(self):
        """
        The right dict must be returned when a long structure sequence is
        passed.
        """
        structure = '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---'
        landmark = PredictedStructure()
        result = landmark.getStructureOffsets(structure)
        self.assertEqual({'-': [[0, 0], [28, 28], [31, 31], [36, 38], [42, 42],
                                [49, 49], [52, 54]],
                          'H': [[1, 26], [45, 48]],
                          'S': [[34, 35], [50, 51]],
                          'T': [[27, 27], [29, 30], [32, 33], [39, 41],
                                [43, 44]]}, result)


class TestPredictedStructureFind(TestCase):
    """
    Tests for the PredictedStructure.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testIdenticalStructureRead(self):
        """
        A read that only consists of one type of structure must lead to the
        correct landmark.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGGGGGGG')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual([Landmark('G', 'ST', 0, 9)], list(result))

    def testFirstItemDifferentRead(self):
        """
        If the first item in the structure sequence is different from the rest,
        the right landmarks must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HGGGGGGGG')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual(set([Landmark('G', 'ST', 1, 8),
                              Landmark('H', 'ST', 0, 1)]), set(result))

    def testLastItemDifferentRead(self):
        """
        If the last item in the structure sequence is different from the rest,
        the right landmarks must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGGGGGGH')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual(set([Landmark('H', 'ST', 8, 1),
                              Landmark('G', 'ST', 0, 8)]), set(result))

    def testLongStructureSequenceNotDefault(self):
        """
        The right landmarks must be returned when a long structure sequence is
        passed in, but most of the structures are not specified in
        structureNames.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual(set([Landmark('H', 'ST', 1, 26),
                              Landmark('H', 'ST', 45, 4)]), set(result))

    def testLongStructureSequenceDefault(self):
        """
        The right landmarks must be returned when a long structure sequence is
        passed in, and most of the structures are  specified in structureNames.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-GG-TTSS---III-IIHHHH-EEEEE')
        landmark = PredictedStructure()
        result = landmark.find(read)
        self.assertEqual(set([Landmark('H', 'ST', 1, 26),
                              Landmark('H', 'ST', 45, 4),
                              Landmark('G', 'ST', 29, 2),
                              Landmark('I', 'ST', 39, 3),
                              Landmark('I', 'ST', 43, 2),
                              Landmark('E', 'ST', 50, 5)]), set(result))
