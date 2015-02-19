from unittest import TestCase

from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn)
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs)


from light.performance.query import queryDatabase


class _TestPolymerase(object):
    """
    Test look-up of polymerase sequences.
    """

    landmarkClasses = None  # Must be set in a subclass.
    trigClasses = None  # Must be set in a subclass.
    limitPerLandmark = None
    maxDistance = None
    minDistance = None

    def testFind(self):
        """
        Look up various polymerase sequences.
        """
        self.details = queryDatabase(
            'performance/database/polymerase-db.fasta',
            'performance/read/polymerase-queries.fasta',
            self.landmarkClasses, self.trigClasses,
            self.limitPerLandmark, self.maxDistance, self.minDistance)


class TestAlphaHelix(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix landmark finder.
    """
    landmarkClasses = [AlphaHelix]
    trigClasses = []


class TestAlphaHelix_3_10(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix_3_10 landmark finder.
    """
    landmarkClasses = [AlphaHelix_3_10]
    trigClasses = []


class TestAlphaHelix_pi(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix_pi landmark finder.
    """
    landmarkClasses = [AlphaHelix_pi]
    trigClasses = []


class TestBetaStrand(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaStrand landmark finder.
    """
    landmarkClasses = [BetaStrand]
    trigClasses = []


class TestBetaTurn(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaTurn landmark finder.
    """
    landmarkClasses = [BetaTurn]
    trigClasses = []


class TestAminoAcids(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AminoAcids trig point finder.
    """
    landmarkClasses = []
    trigClasses = [AminoAcids]


class TestPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Peaks trig point finder.
    """
    landmarkClasses = []
    trigClasses = [Peaks]


class TestTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Troughs trig point finder.
    """
    landmarkClasses = []
    trigClasses = [Troughs]


class TestIndividualPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an IndividualPeaks trig point finder.
    """
    landmarkClasses = []
    trigClasses = [IndividualPeaks]


class TestIndividualTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an IndividualTroughs trig point
    finder.
    """
    landmarkClasses = []
    trigClasses = [IndividualTroughs]
