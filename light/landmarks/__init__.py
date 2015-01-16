from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi
from light.landmarks.amino_acids import AminoAcids
from light.landmarks.beta_strand import BetaStrand
from light.landmarks.beta_turn import BetaTurn
from light.landmarks.gor4_alpha_helix import GOR4AlphaHelix
from light.landmarks.gor4_beta_strand import GOR4BetaStrand
from light.landmarks.prosite import Prosite

ALL_LANDMARK_FINDER_CLASSES = {
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids, BetaStrand,
    BetaTurn, GOR4AlphaHelix, GOR4BetaStrand, Prosite}

DEFAULT_LANDMARK_FINDER_CLASSES = {
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids, BetaStrand,
    BetaTurn, Prosite}


def findLandmark(name):
    """
    A function to find a landmark finder.

    @param name: The C{str} name of the landmark finder class to find.
    @return: The found class, or C{None}.
    """

    for klass in ALL_LANDMARK_FINDER_CLASSES:
        if name == klass.NAME:
            return klass


# Default exports for 'from light.landmarks import *'
__all__ = ['findLandmark', 'ALL_LANDMARK_FINDER_CLASSES',
           'DEFAULT_LANDMARK_FINDER_CLASSES']
