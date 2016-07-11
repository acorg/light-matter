from light.features import Landmark
from light.finder import Finder
from light.utils import stringSpans


def combineHelices(structure):
    """
    A function which returns a string where all amino acids assigned 'H', 'G'
    or 'I' are given the same letter.

    @param structure: A C{str} of structure information.
    """
    combinedStructure = ''
    for aa in structure:
        if aa in {'H', 'G', 'I'}:
            combinedStructure += 'K'
        else:
            combinedStructure += 'C'

    return combinedStructure


class PDB_CombinedAlphaHelix(Finder):
    """
    A class for finding PDB alpha helices based on predicted secondary
    structures. All three alpha helix types (AlphaHelix, AlphaHelix_3_10 and
    AlphaHelix_pi) are considered together.

    Predicted secondary structure sequences can be downloaded from
    http://www.rcsb.org/pdb/files/ss.txt. The PDB predictions are made using
    the DSSP algorithm.
    """
    NAME = 'PDB Combined AlphaHelix'
    SYMBOL = 'PDB-C-A'
    STRUCTURE_LETTER = 'K'

    def find(self, read):
        """
        Find landmarks based on known secondary structures.

        @param read: An instance of C{dark.reads.SSAARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        combinedStructure = combineHelices(read.structure)
        for letter, start, end in stringSpans(combinedStructure):
            if letter == self.STRUCTURE_LETTER:
                yield Landmark(self.NAME, self.SYMBOL, start, end - start)
