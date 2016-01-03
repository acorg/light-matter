from collections import defaultdict

from light.features import Landmark, Finder


class PredictedStructure(Finder):
    """
    A class for computing statistics based on predicted secondary structures.
    The predicted secondary structure sequences were downloaded from
    http://www.rcsb.org/pdb/files/ss.txt on the 11/11/2015. The predictions are
    made using the DSSP algorithm.
    This finder is intended to test the performance of our algorithm if perfect
    secondary structure information is available to be used as landmarks.
    """
    NAME = 'PredictedStructure'
    SYMBOL = 'ST'

    def getStructureOffsets(self, ssSequence):
        """
        Takes a structure sequence, and returns a dictionary that contains the
        offsets for each structure.

        @param ssSequence: A C{str} structure sequence.
        """
        result = defaultdict(list)

        previous = None
        start = 0
        for i, item in enumerate(ssSequence):
            # item is the last item of the sequence
            try:
                ssSequence[i + 1]
            except IndexError:
                if item == previous:
                    result[item].append([start, i])
                else:
                    result[item].append([i, i])
                break
            if item == previous and ssSequence[i + 1] != item:
                result[item].append([start, i])
            # item is the only item of the sequence
            elif item != previous and ssSequence[i + 1] != item:
                previous = item
                result[item].append([i, i])
            # item is the first one in the sequence:
            elif item != previous and ssSequence[i + 1] == item:
                previous = item
                start = i
        return result

    def find(self, read, structureNames=None):
        """
        A function that checks if and where an alpha helix in a sequence
        occurs.

        @param read: An instance of C{light.performance.overlap.SSAARead}.
        @param structureNames: A C{list} of structure names that should be used
            as finders.
        """
        structureNames = structureNames or ['H', 'G', 'I', 'E']

        structureOffsets = self.getStructureOffsets(read.structure)

        for structureName, offsets in structureOffsets.items():
            if structureName in structureNames:
                for offsetPair in offsets:
                    start = offsetPair[0]
                    end = offsetPair[1]
                    yield Landmark(structureName, self.SYMBOL, start,
                                   end - start + 1)
