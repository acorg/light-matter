from .pdb import PDB_Finder


class PDB_AlphaHelix(PDB_Finder):
    """
    Find PDB Alpha Helix landmarks.

    See C{light.landmarks.pdb.PDB_Finder} for details.
    """
    NAME = 'PDB AlphaHelix'
    SYMBOL = 'PDB-A'
    STRUCTURE_LETTER = 'H'
