from .pdb import PDB_Finder


class PDB_AlphaHelix_pi(PDB_Finder):
    """
    Find PDB Alpha Helix pi landmarks.

    See C{light.landmarks.pdb.PDB_Finder} for details.
    """
    NAME = 'PDB AlphaHelix_pi'
    SYMBOL = 'PDB-Api'
    STRUCTURE_LETTER = 'I'
