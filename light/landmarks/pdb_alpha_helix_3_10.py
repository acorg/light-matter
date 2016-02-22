from .pdb import PDB_Finder


class PDB_AlphaHelix_3_10(PDB_Finder):
    """
    Find PDB Alpha Helix 3, 10 landmarks.

    See C{light.landmarks.pdb.PDB_Finder} for details.
    """
    NAME = 'PDB AlphaHelix_3_10'
    SYMBOL = 'PDB-A310'
    STRUCTURE_LETTER = 'G'
