from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi
from light.landmarks.amino_acids import AminoAcids
from light.landmarks.cluster_alpha_helix import ClusterAlphaHelix
from light.landmarks.beta_strand import BetaStrand
from light.landmarks.beta_turn import BetaTurn
from light.landmarks.gor4_alpha_helix import GOR4AlphaHelix
from light.landmarks.gor4_beta_strand import GOR4BetaStrand
from light.landmarks.gor4_coil import GOR4Coil
from light.landmarks.pdb_alpha_helix import PDB_AlphaHelix
from light.landmarks.pdb_alpha_helix_3_10 import PDB_AlphaHelix_3_10
from light.landmarks.pdb_alpha_helix_pi import PDB_AlphaHelix_pi
from light.landmarks.pdb_extended_strand import PDB_ExtendedStrand
from light.landmarks.prosite import Prosite
from light.landmarks.random import RandomLandmark
from light.landmarks.th_alpha_helix import THAlphaHelix


ALL_LANDMARK_CLASSES = [
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids, BetaStrand,
    BetaTurn, GOR4AlphaHelix, GOR4BetaStrand, GOR4Coil, Prosite,
    THAlphaHelix, ClusterAlphaHelix]

DEV_LANDMARK_CLASSES = [
    PDB_AlphaHelix, PDB_AlphaHelix_3_10, PDB_AlphaHelix_pi,
    PDB_ExtendedStrand, RandomLandmark]

ALL_LANDMARK_CLASSES_INCLUDING_DEV = (ALL_LANDMARK_CLASSES +
                                      DEV_LANDMARK_CLASSES)

DEFAULT_LANDMARK_CLASSES = [
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids, BetaStrand,
    BetaTurn, GOR4AlphaHelix, GOR4BetaStrand, Prosite]

_HASHKEY_TO_NAME = dict((cls.SYMBOL, cls.NAME) for cls in ALL_LANDMARK_CLASSES)


def findLandmark(name):
    """
    A function to find a landmark finder by name.

    @param name: The C{str} name of the landmark finder class to find.
    @return: The found class, or C{None}.
    """

    for klass in ALL_LANDMARK_CLASSES:
        if name == klass.NAME:
            return klass


def findLandmarks(names):
    """
    A function to find multiple landmark finders by name.

    @param names: A C{list} of C{str} name of the landmark finder classes to
        find, or C{None} (equivalent to an empty list of names)..
    @raise ValueError: If any name cannot be found.
    @return: A C{list} of classes that were successfully found.
    """
    found = []
    unknown = []
    if names:
        for name in names:
            klass = findLandmark(name)
            if klass:
                found.append(klass)
            else:
                unknown.append(name)
    if unknown:
        raise ValueError('Unknown landmark finder%s: %s.' % (
            's' if len(unknown) > 1 else '', ', '.join(unknown)))
    else:
        return found


def landmarkNameFromHashkey(hashkey):
    """
    Get the name of a landmark class given part of a hashkey that was
    created for it (in features.Landmark.hashkey).

    @param hashkey: A C{str} hashkey for the class.
    @return: A C{str} landmark class name, or C{None} if the hashkey
        cannot be found.
    """
    # Look for the entire hashkey in _HASHKEY_TO_NAME and if that's not
    # present, shorten it one character at a time from the right. This
    # allows us to find the correct finder class for hashkeys like PS00342
    # (the Prosite class, whose symbol is PS).
    while hashkey:
        try:
            return _HASHKEY_TO_NAME[hashkey]
        except KeyError:
            hashkey = hashkey[:-1]


# Default exports for 'from light.landmarks import *'
__all__ = ['findLandmark', 'findLandmarks', 'landmarkNameFromHashkey',
           'ALL_LANDMARK_CLASSES', 'DEFAULT_LANDMARK_CLASSES']
