from .pdb import PDB_Finder


class PDB_ExtendedStrand(PDB_Finder):
    """
    Find PDB extended strand landmarks.

    See C{light.landmarks.pdb.PDB_Finder} for details.
    """
    NAME = 'PDB ExtendedStrand'
    SYMBOL = 'PDB-ES'
    STRUCTURE_LETTER = 'E'
