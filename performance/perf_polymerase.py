from unittest import TestCase

from light.database import Database
from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn)
from light.parameters import DatabaseParameters
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs)


from light.performance.query import queryDatabase


class _TestPolymerase(object):
    """
    Test look-up of polymerase sequences.
    """

    LANDMARKS = None  # Must be set in a subclass.
    TRIG_POINTS = None  # Must be set in a subclass.
    LIMIT_PER_LANDMARK = None
    MAX_DISTANCE = None
    MIN_DISTANCE = None

    def testFind(self):
        """
        Look up various polymerase sequences.
        """
        dbParams = DatabaseParameters(landmarks=self.LANDMARKS,
                                      trigPoints=self.TRIG_POINTS,
                                      limitPerLandmark=self.LIMIT_PER_LANDMARK,
                                      maxDistance=self.MAX_DISTANCE,
                                      minDistance=self.MIN_DISTANCE)
        database = Database(dbParams=dbParams)
        self.details = queryDatabase(
            'performance/database/polymerase-db.fasta',
            'performance/read/polymerase-queries.fasta',
            database)


class TestAlphaHelix(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix landmark finder.
    """
    LANDMARKS = [AlphaHelix]
    TRIG_POINTS = []


class TestAlphaHelix_3_10(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix_3_10 landmark finder.
    """
    LANDMARKS = [AlphaHelix_3_10]
    TRIG_POINTS = []


class TestAlphaHelix_pi(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AlphaHelix_pi landmark finder.
    """
    LANDMARKS = [AlphaHelix_pi]
    TRIG_POINTS = []


class TestBetaStrand(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaStrand landmark finder.
    """
    LANDMARKS = [BetaStrand]
    TRIG_POINTS = []


class TestBetaTurn(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a BetaTurn landmark finder.
    """
    LANDMARKS = [BetaTurn]
    TRIG_POINTS = []


class TestAminoAcids(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an AminoAcids trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [AminoAcids]


class TestPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Peaks trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [Peaks]


class TestTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just a Troughs trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [Troughs]


class TestIndividualPeaks(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an IndividualPeaks trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [IndividualPeaks]


class TestIndividualTroughs(_TestPolymerase, TestCase):
    """
    Test looking up polymerases with just an IndividualTroughs trig point
    finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [IndividualTroughs]
