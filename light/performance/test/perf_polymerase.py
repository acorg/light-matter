from os import mkdir
from os.path import dirname, exists, isdir, join
from unittest import TestCase
import matplotlib.pyplot as plt
from scipy import stats

from dark.fasta_ss import SSFastaReads

import light
from light.database import Database
from light.performance.affinity import affinityMatrix
from light.performance.polymerase import Z_SCORES, BIT_SCORES
from light.performance import testArgs

# Create a singleton affinity matrix of lm scores for all polymerase
# sequences.
_QUERIES = list(SSFastaReads(
    join(dirname(light.__file__),
         'performance', 'data', 'polymerase.fasta')))
_AFFINITY = affinityMatrix(_QUERIES, database=Database(testArgs.dbParams),
                           findParams=testArgs.findParams, returnDict=True)


def makeOutputDir(scoreTypeX, scoreTypeY):
    """
    Create an output directory if it doesn't already exist.

    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @return: The C{str} path to the output directory.
    """
    FILESYSTEM_NAME = {
        'Bit score': 'bit-score',
        'Z score': 'z-score',
        'Light matter score': 'lm-score',
    }

    assert scoreTypeX in FILESYSTEM_NAME and scoreTypeY in FILESYSTEM_NAME

    def makeDir(path):
        """
        Check to see if a path exists and whether it's a directory. Create
        it (as a dir) if it's non-existent.

        @raise RuntimeError: if C{path} exists but is not a directory.
        @param path: The C{str} path for the directory.
        """
        if exists(path):
            if not isdir(path):
                raise RuntimeError('Output path %r already exists but is not '
                                   'a directory.' % path)
        else:
            mkdir(path)

    outputDir = join(testArgs.outputDir, 'polymerase')

    subDir = join(outputDir,
                  '%s-vs-%s' % (FILESYSTEM_NAME[scoreTypeX],
                                FILESYSTEM_NAME[scoreTypeY]))
    makeDir(outputDir)
    makeDir(subDir)

    return subDir


def plot(x, y, readId, scoreTypeX, scoreTypeY):
    """
    Make a scatterplot of the test results.

    @param x: a C{list} of C{float} x coordinates.
    @param y: a C{list} of C{float} y coordinates.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    """
    MAX_Z_SCORE = 60.0

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
                 (readId, rValue, se, slope, pValue))
    ax.set_ylabel(scoreTypeY)
    ax.set_xlabel(scoreTypeX)

    # Light matter scores are always 0.0 to 1.0 so we can explicitly set
    # axis limits for these.
    if scoreTypeX == 'Light matter score':
        ax.set_xlim(0.0, 1.0)
    if scoreTypeY == 'Light matter score':
        ax.set_ylim(0.0, 1.0)

    # Z scores are always 0.0 to MAX_Z_SCORE so we can explicitly set
    # axis limits for these.
    if scoreTypeX == 'Z score':
        ax.set_xlim(0.0, MAX_Z_SCORE)
    if scoreTypeY == 'Z score':
        ax.set_ylim(0.0, MAX_Z_SCORE)

    # Axes.
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    fig.savefig(join(makeOutputDir(scoreTypeX, scoreTypeY),
                     '%s.png' % readId))
    plt.close()


class TestZScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        for queryId in Z_SCORES:
            zScores = []
            lmScores = []
            for subjectId in Z_SCORES:
                if queryId != subjectId:
                    lmScores.append(_AFFINITY[queryId][subjectId])
                    zScores.append(Z_SCORES[queryId][subjectId])
            plot(lmScores, zScores, queryId, 'Light matter score', 'Z score')


class TestBitScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        for queryId in BIT_SCORES:
            zScores = []
            lmScores = []
            for subjectId in BIT_SCORES:
                if queryId != subjectId:
                    lmScores.append(_AFFINITY[queryId][subjectId])
                    zScores.append(BIT_SCORES[queryId][subjectId])
            plot(lmScores, zScores, queryId, 'Light matter score', 'Bit score')


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between Z scores and blast bit scores.
        """
        for queryId in BIT_SCORES:
            zScores = []
            bitScores = []
            for subjectId in BIT_SCORES:
                if queryId != subjectId:
                    bitScores.append(BIT_SCORES[queryId][subjectId])
                    zScores.append(Z_SCORES[queryId][subjectId])
            plot(bitScores, zScores, queryId, 'Bit score', 'Z score')
