from light.features import Landmark
from light.finder import Finder
from light.utils import stringSpans


class PDB_Finder(Finder):
    """
    A class for finding PDB landmarks based on predicted secondary structures.

    This class is not intended to be used directly, but to be subclassed.

    Predicted secondary structure sequences can be
    downloaded from http://www.rcsb.org/pdb/files/ss.txt. The PDB predictions
    are made using the DSSP algorithm.  Subclassed of this finder are used to
    test the performance of our algorithm if perfect secondary structure
    information is available from PDB to use as landmarks.
    """
    # The following three variables must be set in subclasses.
    NAME = SYMBOL = STRUCTURE_LETTER = None

    def find(self, read):
        """
        Find landmarks based on known secondary structures.

        @param read: An instance of C{dark.reads.SSAARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for letter, start, end in stringSpans(read.structure):
            if letter == self.STRUCTURE_LETTER:
                yield Landmark(self.NAME, self.SYMBOL, start, end - start)