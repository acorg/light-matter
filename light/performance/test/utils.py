from os import mkdir
from os.path import exists, isdir, join
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from light.performance import testArgs

# Keep pyflakes quiet by pretending to use Axes3D.
_ = Axes3D


def makeDir(path):
    """
    Check to see if a path exists and whether it's a directory. Create it (as
    a dir) if it's non-existent.

    @raise RuntimeError: if C{path} exists but is not a directory.
    @param path: The C{str} path for the directory.
    """
    if exists(path):
        if not isdir(path):
            raise RuntimeError('Output path %r already exists but is not a '
                               'directory.' % path)
    else:
        mkdir(path)


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


def plot3D(x, y, z, readId, scoreTypeX, scoreTypeY, scoreTypeZ,
           interactive=False):
    """
    Make a 3D plot of the test results.

    @param x: a C{list} of C{float} X axis bit score values.
    @param y: a C{list} of C{float} Y axis Z score values.
    @param z: a C{list} of C{float} Z axis light matter score values.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param scoreTypeZ: A C{str} Z-axis title indicating the type of score.
    @param interactive: If C{True} use plt.show() to display interactive plots
        that the user will need to manually dismiss.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    alpha = 0.1
    dotSize = 40

    # The initial view onto the 3D plot.
    cameraDegrees = 11
    azimuth = -103
    ax.view_init(cameraDegrees, azimuth)

    # All the Z scores we've observed so far are less than 65. We have a
    # fixed upper limit here so all graphs produced will have the same Y
    # axis upper limit.
    zScoreLimit = 65.0
    assert max(y) <= zScoreLimit

    # Values less than these cutoffs are considered bad, values bigger are
    # good. Planes will be drawn to separate each axis into bad/good
    # points (assuming there are any good points in a given dimension).
    bitScoreCutoff = 50.0
    zScoreCutoff = 20.0
    lmScoreCutoff = 0.5

    # cutoffStringToColor converts a binary string of length 3 to a color.
    # The positions in the string represent bad/good (0/1) values for bit
    # score, Z score, light matter score according to the bad/good cutoff
    # values above. From best to worst: yellow, green, black, blue, red.
    # Blue isn't bad, it just indicates that two sequences are similar but
    # that they have no common structure (they may have no structure at
    # all).
    cutoffStringToColor = {
        '000': 'green',   # OK - no conflict.
        '001': 'red',     # Really bad - lm disagrees with both other scores.
        '010': 'black',   # Quite bad - Lm disagrees with Z score.
        '011': 'yellow',  # Good - bit score low, structure scores both high.
        '100': 'blue',    # Weird - bit score is the only one that's high.
        '101': 'black',   # Quite bad - Lm disagrees with Z score.
        '110': 'red',     # Really bad - lm disagrees with both other scores.
        '111': 'green',   # OK - no conflict.
    }

    # Assign each x, y, z triple a color.
    colors = []
    for bitScore, zScore, lmScore in zip(x, y, z):
        key = (('1' if bitScore > bitScoreCutoff else '0') +
               ('1' if zScore > zScoreCutoff else '0') +
               ('1' if lmScore > lmScoreCutoff else '0'))
        colors.append(cutoffStringToColor[key])

    ax.scatter(x, y, z, c=colors, s=dotSize)
    ax.set_xlabel(scoreTypeX)
    ax.set_ylabel(scoreTypeY)
    ax.set_zlabel(scoreTypeZ)
    ax.set_title(readId)

    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, zScoreLimit)
    ax.set_zlim(0.0, 1.0)

    # cg is the contour granularity: the number of points in the mesh used
    # for plotting the three plane contours that divide each axis into good
    # & bad regions.
    cg = 60

    # Bit score (x) bad/good plane.
    if max(x) >= bitScoreCutoff:
        Y = np.linspace(0.0, zScoreLimit, cg)
        Z = np.linspace(0.0, 1.0, cg)
        yy, zz = np.meshgrid(Y, Z)
        X = np.array([bitScoreCutoff] * (cg * cg)).reshape(cg, cg)
        ax.contourf(X, yy, zz, colors='blue', alpha=alpha, zdir='x')

    # Z score (y) bad/good plane.
    if max(y) >= zScoreCutoff:
        X = np.linspace(0.0, max(x), cg)
        Z = np.linspace(0.0, 1.0, cg)
        xx, zz = np.meshgrid(X, Z)
        Y = np.array([zScoreCutoff] * (cg * cg)).reshape(cg, cg)
        ax.contourf(xx, Y, zz, colors='green', alpha=alpha, zdir='y')

    # LM score (z) bad/good plane.
    if max(z) >= lmScoreCutoff:
        X = np.linspace(0.0, max(x), cg)
        Y = np.linspace(0.0, zScoreLimit, cg)
        xx, yy = np.meshgrid(X, Y)
        Z = np.array([lmScoreCutoff] * (cg * cg)).reshape(cg, cg)
        ax.contourf(xx, yy, Z, colors='purple', alpha=alpha, zdir='z')

    outputDir = join(testArgs.outputDir, 'polymerase')
    subDir = join(outputDir, '3d')
    makeDir(outputDir)
    makeDir(subDir)

    fig.savefig(join(subDir, '%s.png' % readId))
    if interactive:
        plt.show()
