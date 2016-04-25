from os import mkdir
from os.path import exists, isdir, join
import matplotlib.pyplot as plt
from scipy import stats

from light.performance import testArgs


def makeOutputDir(scoreTypeX, scoreTypeY, dirName):
    """
    Create an output directory if it doesn't already exist.

    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param dirName: A C{str} name of the output directory to be created.

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

    outputDir = join(testArgs.outputDir, dirName)

    subDir = join(outputDir,
                  '%s-vs-%s' % (FILESYSTEM_NAME[scoreTypeX],
                                FILESYSTEM_NAME[scoreTypeY]))
    makeDir(outputDir)
    makeDir(subDir)

    return subDir


def plot(x, y, readId, scoreTypeX, scoreTypeY, dirName):
    """
    Make a scatterplot of the test results.

    @param x: a C{list} of C{float} x coordinates.
    @param y: a C{list} of C{float} y coordinates.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param dirName: A C{str} name of the output directory to be created.
    """
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    slope, intercept, rValue, pValue, se = stats.linregress(x, y)

    # Plot.
    plt.plot(x, y, 'o', markerfacecolor='blue', markeredgecolor='blue')
    plt.plot([0, max(x)], [intercept, slope * max(x) + intercept], '-',
             color='green' if slope >= 0 else 'red')

    # Labels.
    ax.set_title('Read: %s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' %
                 (readId, rValue, se, slope, pValue))
    ax.set_ylabel(scoreTypeY)
    ax.set_xlabel(scoreTypeX)

    # Light matter scores are always <= 1.0.
    if scoreTypeX == 'Light matter score':
        ax.set_xlim(right=1.0)
    if scoreTypeY == 'Light matter score':
        ax.set_ylim(top=1.0)

    # Z scores are always <= 60.0 (with sanity check).
    if scoreTypeX == 'Z score':
        assert max(x) <= 65.0
        ax.set_xlim(right=65.0)
    if scoreTypeY == 'Z score':
        assert max(y) <= 65.0
        ax.set_ylim(top=65.0)

    # No scores can be negative. Explicitly set the lower limits on both
    # axes to zero. This stops regression lines from causing an axis to
    # display useless areas with negative ticks and no data points.
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)

    # Axes.
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    fig.savefig(join(makeOutputDir(scoreTypeX, scoreTypeY, dirName),
                     '%s.png' % readId))
    plt.close()
