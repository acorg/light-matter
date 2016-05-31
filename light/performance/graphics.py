"""
None of this code is tested.
"""

from os.path import join
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import plotly
import plotly.graph_objs as go

from light.performance.utils import pythonNameToPdbName

# Keep pyflakes quiet by pretending to use Axes3D.
_ = Axes3D


def plot(x, y, readId, scoreTypeX, scoreTypeY, dirName):
    """
    Make a scatterplot of the test results.

    @param x: a C{list} of C{float} x coordinates.
    @param y: a C{list} of C{float} y coordinates.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param dirName: A C{str} name of the output directory in which to store
        the plot image. The image will be saved to dirName + readId + '.png'
    @raises AssertionError: If the length of C{x} is not the same as the
        length of C{y}.
    """
    assert len(x) == len(y)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    if len(x) > 1:
        slope, intercept, rValue, pValue, se = stats.linregress(x, y)

        # Plot.
        plt.plot(x, y, 'o', markerfacecolor='blue', markeredgecolor='white')
        plt.plot([0, max(x)], [intercept, slope * max(x) + intercept], '-',
                 color='green' if slope >= 0 else 'red')

        # Labels.
        ax.set_title('Read: %s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' %
                     (pythonNameToPdbName(readId), rValue, se, slope, pValue))
    else:
        ax.set_title('No (or not enough) x,y data given for read %s' %
                     pythonNameToPdbName(readId))

    ax.set_ylabel(scoreTypeY)
    ax.set_xlabel(scoreTypeX)

    # Light matter scores are always <= 1.0.
    if scoreTypeX == 'Light matter score':
        ax.set_xlim(right=1.0)
    if scoreTypeY == 'Light matter score':
        ax.set_ylim(top=1.0)

    # Z scores are always <= 60.0 (with sanity check).
    if scoreTypeX == 'Z score':
        if x:
            assert max(x) <= 65.0
        ax.set_xlim(right=65.0)
    if scoreTypeY == 'Z score':
        if y:
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

    fig.savefig(join(dirName, '%s.png' % readId))
    plt.close()


def plot3D(x, y, z, readId, scoreTypeX, scoreTypeY, scoreTypeZ,
           dirName, interactive=False):
    """
    Make a 3D plot of the test results.

    @param x: a C{list} of C{float} X axis bit score values.
    @param y: a C{list} of C{float} Y axis Z score values.
    @param z: a C{list} of C{float} Z axis light matter score values.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param scoreTypeZ: A C{str} Z-axis title indicating the type of score.
    @param dirName: A C{str} name of the output directory in which to store
        the plot image. The image will be saved to dirName + readId + '.png'
    @param interactive: If C{True} use plt.show() to display interactive plots
        that the user will need to manually dismiss.
    @raises AssertionError: If the length of C{x} is not the same as the
        length of C{y} and the length of C{z}.
    """
    assert len(x) == len(y) == len(z)

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
    if y:
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

    if x:
        # Assign each x, y, z triple a color.
        colors = []
        for bitScore, zScore, lmScore in zip(x, y, z):
            key = (('1' if bitScore > bitScoreCutoff else '0') +
                   ('1' if zScore > zScoreCutoff else '0') +
                   ('1' if lmScore > lmScoreCutoff else '0'))
            colors.append(cutoffStringToColor[key])

        ax.scatter(x, y, z, c=colors, s=dotSize)
        ax.set_title(pythonNameToPdbName(readId))

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
    else:
        ax.set_title('No x,y,z data given for read %s' %
                     pythonNameToPdbName(readId))

    ax.set_xlabel(scoreTypeX)
    ax.set_ylabel(scoreTypeY)
    ax.set_zlabel(scoreTypeZ)

    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, zScoreLimit)
    ax.set_zlim(0.0, 1.0)

    fig.savefig(join(dirName, '%s.png' % readId))
    if interactive:
        plt.show()
    plt.close()


def plot3DPlotly(x, y, z, readId, scoreTypeX, scoreTypeY, scoreTypeZ,
                 dirName, interactive=False):
    """
    Make a 3D plot of the test results using Plotly.

    @param x: a C{list} of C{float} X axis bit score values.
    @param y: a C{list} of C{float} Y axis Z score values.
    @param z: a C{list} of C{float} Z axis light matter score values.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param scoreTypeZ: A C{str} Z-axis title indicating the type of score.
    @param dirName: A C{str} name of the output directory in which to store
        the plot image. The image will be saved to dirName + readId + '.png'
    @param interactive: If C{True} use plt.show() to display interactive plots
        that the user will need to manually dismiss.
    @raises AssertionError: If the length of C{x} is not the same as the
        length of C{y} and the length of C{z}.
    """
    assert len(x) == len(y) == len(z)

    trace1 = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=12,
            line=dict(
                color='rgba(217, 217, 217, 0.14)',
                width=0.5
            ),
            opacity=0.8
        )
    )

    data = [trace1]
    layout = go.Layout(
        title=pythonNameToPdbName(readId),
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        xaxis=dict(
            title=scoreTypeX,
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        ),
        yaxis=dict(
            title=scoreTypeY,
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    filename = join(dirName, '%s.html' % readId)
    plotly.offline.plot(fig, show_link=False, filename=filename,
                        auto_open=interactive)
