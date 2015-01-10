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

    landmarkFinderClasses = None  # Must be set in a subclass.
    trigFinderClasses = None  # Must be set in a subclass.
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
            self.landmarkFinderClasses, self.trigFinderClasses,
            self.limitPerLandmark, self.maxDistance, self.minDistance)


class TestAlphaHelix(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a AlphaHelix landmark finder.
    """
    landmarkFinderClasses = [AlphaHelix]
    trigFinderClasses = []


class TestAlphaHelix_3_10(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a AlphaHelix_3_10 landmark finder.
    """
    landmarkFinderClasses = [AlphaHelix_3_10]
    trigFinderClasses = []


class TestAlphaHelix_pi(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a AlphaHelix_pi landmark finder.
    """
    landmarkFinderClasses = [AlphaHelix_pi]
    trigFinderClasses = []


class TestBetaStrand(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaStrand landmark finder.
    """
    landmarkFinderClasses = [BetaStrand]
    trigFinderClasses = []


class TestBetaTurn(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaTurn landmark finder.
    """
    landmarkFinderClasses = [BetaTurn]
    trigFinderClasses = []


class TestAminoAcids(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a AminoAcids trig point finder.
    """
    landmarkFinderClasses = []
    trigFinderClasses = [AminoAcids]


class TestPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Peaks trig point finder.
    """
    landmarkFinderClasses = []
    trigFinderClasses = [Peaks]


class TestTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Troughs trig point finder.
    """
    landmarkFinderClasses = []
    trigFinderClasses = [Troughs]


class TestIndividualPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a IndividualPeaks trig point finder.
    """
    landmarkFinderClasses = []
    trigFinderClasses = [IndividualPeaks]


class TestIndividualTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a IndividualTroughs trig point
    finder.
    """
    landmarkFinderClasses = []
    trigFinderClasses = [IndividualTroughs]
