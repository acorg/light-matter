import time
from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn)
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs)


class TestDatabaseCreation(TestCase):
    """
    Test database creation performance.
    """

    def testCreation(self):
        """
        How long does it take to construct a database with no subjects, and no
        landmark or trig point finders?
        """
        Database([], [])


class _SubjectAddMixin(object):
    """
    Test database subject adding performance.
    """

    # The number of copies of a known feature that should be present in the
    # sequence in the read that is added to the database. See our subclasses
    # for usage.
    FEATURE_COUNT = 50

    # The number of subjects to add to the database.
    SUBJECT_COUNT = 50

    def testAddSubjects(self):
        """
        How long does it take to add SUBJECT_COUNT subjects to a database.
        """
        read = AARead('id', self.FEATURE_SEQUENCE * self.FEATURE_COUNT)
        database = Database(self.LANDMARKS, self.TRIG_POINTS)
        startTime = time.time()
        for _ in xrange(self.SUBJECT_COUNT):
            database.addSubject(read)
        elapsed = time.time() - startTime
        self.details = elapsed


class SubjectAddNoFinders(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has no landmark or
    trig point finders.
    """
    LANDMARKS = []
    TRIG_POINTS = []
    # The sequence below be random. It's not going to be examined by
    # anything because there are no landmark or trig point finders.
    FEATURE_SEQUENCE = 'MMMMMMMMMM'


class SubjectAddBetaStrand(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has a BetaStrand
    landmark finder.
    """
    LANDMARKS = [BetaStrand]
    TRIG_POINTS = []
    SEP = 'N'
    FEATURE_SEQUENCE = BetaStrand.BETA_STRAND_AAs + SEP

    def testSeparatorIsntABetaStrandResidue(self):
        """
        The residue we use to separate Beta strands in the test sequence
        must not itself be a Beta strand residue.
        """
        self.assertNotIn(self.SEP, BetaStrand.BETA_STRAND_AAs)
        self.ignore = True


class SubjectAddBetaTurn(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has a BetaTurn
    landmark finder.
    """
    LANDMARKS = [BetaTurn]
    TRIG_POINTS = []
    FEATURE_SEQUENCE = 'NPNWA'


class SubjectAddAlphaHelix(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an AlphaHelix
    landmark finder.
    """
    LANDMARKS = [AlphaHelix]
    TRIG_POINTS = []
    FEATURE_SEQUENCE = 'FRRRFRRRFA'


class SubjectAddAlphaHelix_3_10(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    AlphaHelix_3_10 landmark finder.
    """
    LANDMARKS = [AlphaHelix_3_10]
    TRIG_POINTS = []
    FEATURE_SEQUENCE = 'FRRFRRFA'


class SubjectAddAlphaHelix_pi(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    AlphaHelix_pi landmark finder.
    """
    LANDMARKS = [AlphaHelix_pi]
    TRIG_POINTS = []
    FEATURE_SEQUENCE = 'FRRRRFRRRRFA'


class SubjectAddAminoAcids(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    AminoAcids trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [AminoAcids]
    FEATURE_SEQUENCE = 'CMMMMMMM'


class SubjectAddIndividualPeaks(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    IndividualPeaks trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [IndividualPeaks]
    FEATURE_SEQUENCE = 'ATAAAA'


class SubjectAddIndividualTroughs(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    IndividualTroughs trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [IndividualTroughs]
    FEATURE_SEQUENCE = 'TMTTTT'


class SubjectAddPeaks(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    Peaks trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [Peaks]
    FEATURE_SEQUENCE = 'ASAAAA'


class SubjectAddTroughs(_SubjectAddMixin, TestCase):
    """
    A class to test adding subjects to a database that just has an
    Troughs trig point finder.
    """
    LANDMARKS = []
    TRIG_POINTS = [Troughs]
    FEATURE_SEQUENCE = 'AVAAAA'
