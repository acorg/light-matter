from unittest import TestCase
import matplotlib.pyplot as plt
from scipy import stats
from os.path import join

from light.database import Database
from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn)
from light.parameters import DatabaseParameters
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs)

from light.performance.affinity import affinityMatrix
from light.performance.query import queryDatabase
from light.performance.polymerase import Z_SCORES
from light.performance import testArgs


def plot(x, y, read, scoreType):
    """
    Make a scatterplot of the test results.

    @param x: a C{list} of the x coordinates (light matter score).
    @param y: a C{list} of the y coordinates (either z-score or bit score).
    @param scoreType: A C{str} Y-axis title indicating the type of score.
    """
    plt.rcParams['font.family'] = 'Helvetica'
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    slope, intercept, rValue, pValue, se = stats.linregress(x, y)

    # Plot.
    plt.plot(x, y, 'o', markerfacecolor='blue', markeredgecolor='blue')
    if slope >= 0:
        col = 'green'
    else:
        col = 'red'
    plt.plot([0, max(x)], [intercept, slope * max(x) + intercept], '-',
             color=col)

    # Labels.
    ax.set_title('Read: %s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' %
                 (read, rValue, se, slope, pValue))
    ax.set_ylabel(scoreType)
    ax.set_xlabel('Light matter score')

    # Axes.
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    fig.savefig(join(testArgs.outputDir, '%s-%s.png' % (read, scoreType)))
    plt.close()


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


class TestZScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        affinity = affinityMatrix(
            'performance/database/polymerase-db.fasta',
            database=Database(testArgs.dbParams),
            findParams=testArgs.findParams, returnDict=True)

        # Prepare results for plotting, and plot.
        for queryId in Z_SCORES:
            zScores = []
            lmScores = []
            for subjectId in Z_SCORES:
                if queryId != subjectId:
                    lmScores.append(affinity[queryId][subjectId])
                    zScores.append(Z_SCORES[queryId][subjectId])
            plot(lmScores, zScores, queryId, 'Z-score')

        self.details = affinity


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
