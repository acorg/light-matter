import matplotlib.pyplot as plt
from os.path import join
from scipy import stats
from unittest import TestCase

from dark.reads import AARead
from light.database import Database
from light.performance import testArgs
from light.performance.affinity import affinityMatrix
from zscores.hla import BIT_SCORES, Z_SCORES


def plot(x, y, scoreType1, scoreType2=None):
    """
    Make a scatterplot of the test results.

    @param x: A C{list} of the x coordinates (light matter score).
    @param y: A C{list} of the y coordinates (either z-score or bit score).
    @param scoreType1: A C{str} Y-axis title indicating the type of score.
    @param scoreType2: A C{str} X-axis title indicating the type of score. To
        be specified if a score other than the lm score is used.
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
    ax.set_title('Read: 2HLA, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' %
                 (rValue, se, slope, pValue))
    ax.set_ylabel(scoreType1)
    if not scoreType2:
        ax.set_xlabel('Light matter score')
        ax.set_xlim(0.0, 1.0)
    else:
        ax.set_xlabel(scoreType2)

    if scoreType1 == 'Z score':
        ax.set_ylim(0.0, MAX_Z_SCORE)

    # Axes.
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    if not scoreType2:
        fig.savefig(join(testArgs.outputDir, 'hla-%s.png' % (scoreType1)))
    else:
        fig.savefig(join(testArgs.outputDir, 'hla-%s-%s.png' % (scoreType1,
                                                                scoreType2)))

    plt.close()


# Get the light matter scores
subject = AARead('2HLA', 'GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWI'
                         'EQEGPEYWDRNTRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQMMYGCDVG'
                         'SDGRFLRGYRQDAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQW'
                         'RAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSF'
                         'YPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWVAVVVPSGQEQRYTCH'
                         'VQHEGLPKPL')

LM_SCORES = affinityMatrix(
    'performance/read/hla-queries.fasta', subjects=subject,
    database=Database(testArgs.dbParams),
    findParams=testArgs.findParams, returnDict=True)


class TestZScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between our scores and Z scores.
        """
        plot(LM_SCORES, Z_SCORES, 'Z score')


class TestBitScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between our scores and blast bit scores.
        """
        plot(LM_SCORES, BIT_SCORES, 'Bit score')


class TestZScoreBitScoreCorrelation(TestCase):

    def testPlots(self):
        """
        Examine the correlation between Z scores and blast bit scores.
        """
        plot(BIT_SCORES, Z_SCORES, 'Bit score', 'Z score')
