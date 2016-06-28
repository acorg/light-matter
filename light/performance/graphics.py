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

from light.performance.evaluate import PdbSubsetStatistics
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


def plot3DPlotly(bitScores, zScores, lmScores, readId, dirName,
                 scoreTypeX='Bit score', scoreTypeY='Z score',
                 scoreTypeZ='Light matter score', interactive=False,
                 labels=''):
    """
    Make a 3D plot of the test results using Plotly.

    @param bitScores: a C{list} of C{float} bit scores for the X axis.
    @param zScores: a C{list} of C{float} Z scores for the Y axis.
    @param lmScores: a C{list} of C{float} light matter scores for the Z axis.
    @param readId: The C{str} id of the read whose values are being plotted.
    @param dirName: A C{str} name of the output directory in which to store
        the plot image. The image will be saved to dirName + readId + '.png'
    @param scoreTypeX: A C{str} X-axis title indicating the type of score.
    @param scoreTypeY: A C{str} Y-axis title indicating the type of score.
    @param scoreTypeZ: A C{str} Z-axis title indicating the type of score.
    @param interactive: If C{True} use plt.show() to display interactive plots
        that the user will need to manually dismiss.
    @param labels: A C{list} of C{str} labels for the points in the plot. If
        a C{str} is passed, that will be used as the label for all points.
    @raises AssertionError: If the length of C{bitScores} is not the same as
        the length of C{zScores} and the length of C{lmScores}.
    """
    assert len(bitScores) == len(zScores) == len(lmScores)

    # Values less than these cutoffs are considered bad, values bigger are
    # good. Planes will be drawn to separate each axis into good/bad
    # points (assuming there are any good points in a given dimension).
    bitScoreCutoff = 50.0
    zScoreCutoff = 20.0
    lmScoreCutoff = 0.5

    # There's no limit on how high a bit score could be and because their
    # range can be so great it's not practical to have a fixed upper limit
    # for the X axis upper limit (otherwise the plots look weird when no
    # bit scores approach that high upper limit). So we set an artificial
    # max possible bit score based on the values we received and the bit
    # score cut-off value used to display the good/bad bit score plane.
    maxPossibleBitScore = max(max(bitScores), (bitScoreCutoff + 5))
    maxPossibleZScore = 65.0
    maxPossibleLmScore = 1.0

    # All the Z scores we've observed so far are less than 65. We set a
    # fixed upper limit above so all graphs produced will have the same Y
    # axis upper limit.
    if zScores:
        assert max(zScores) <= maxPossibleZScore

    minPossibleBitScore = minPossibleZScore = minPossibleLmScore = 0.0

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

    # Assign each (bit score, Z score, light matter score) triple a color.
    colors = []
    for bitScore, zScore, lmScore in zip(bitScores, zScores, lmScores):
        key = (('1' if bitScore > bitScoreCutoff else '0') +
               ('1' if zScore > zScoreCutoff else '0') +
               ('1' if lmScore > lmScoreCutoff else '0'))
        colors.append(cutoffStringToColor[key])

    # Plot the score triples.
    data = [
        go.Scatter3d(
            x=bitScores,
            y=zScores,
            z=lmScores,
            mode='markers',
            name='Scores',
            marker={
                'size': 9,
                'color': colors,
                'opacity': 0.45,
            },
            text=labels,
        )
    ]

    # The alpha value for the good/bad cut-off planes.
    planeAlpha = 0.1

    planeLine = {
        'color': 'black',
        'width': 1,
    }

    # Bit score (x) good/bad plane.
    data.append(
        go.Scatter3d(
            x=[bitScoreCutoff] * 5,
            y=[minPossibleZScore, maxPossibleZScore, maxPossibleZScore,
               minPossibleZScore, minPossibleZScore],
            z=[maxPossibleLmScore, maxPossibleLmScore, minPossibleLmScore,
               minPossibleLmScore, maxPossibleLmScore],
            mode='lines',
            name=scoreTypeX + ' plane',
            surfaceaxis=0,
            surfacecolor='blue',
            opacity=planeAlpha,
            hoverinfo='none',
            line=planeLine))

    # Z score (y) good/bad plane.
    data.append(
        go.Scatter3d(
            x=[minPossibleBitScore, maxPossibleBitScore,
               maxPossibleBitScore, minPossibleBitScore,
               minPossibleBitScore],
            y=[zScoreCutoff] * 5,
            z=[maxPossibleLmScore, maxPossibleLmScore, minPossibleZScore,
               minPossibleZScore, maxPossibleLmScore],
            mode='lines',
            name=scoreTypeY + ' plane',
            surfaceaxis=1,
            surfacecolor='green',
            opacity=planeAlpha,
            hoverinfo='none',
            line=planeLine))

    # LM score (z) good/bad plane.
    data.append(
        go.Scatter3d(
            x=[minPossibleBitScore, maxPossibleBitScore,
               maxPossibleBitScore, minPossibleBitScore,
               minPossibleBitScore],
            y=[maxPossibleZScore, maxPossibleZScore, minPossibleZScore,
               minPossibleZScore, maxPossibleZScore],
            z=[lmScoreCutoff] * 5,
            mode='lines',
            name=scoreTypeZ + ' plane',
            surfaceaxis=2,
            surfacecolor='purple',
            opacity=planeAlpha,
            hoverinfo='none',
            line=planeLine))

    axisFont = {
        'size': 16,
    }

    layout = go.Layout(
        title=pythonNameToPdbName(readId),
        margin={
            'l': 0,
            'r': 0,
            'b': 0,
            't': 50,
            'pad': 5,
        },
        scene={
            'xaxis': {
                'range': [minPossibleBitScore, maxPossibleBitScore],
                'title': scoreTypeX,
                'titlefont': axisFont,
            },
            'yaxis': {
                'range': [minPossibleZScore, maxPossibleZScore],
                'title': scoreTypeY,
                'titlefont': axisFont,
            },
            'zaxis': {
                'range': [minPossibleLmScore, maxPossibleLmScore],
                'title': scoreTypeZ,
                'titlefont': axisFont,
            },
        },
    )

    fig = go.Figure(data=data, layout=layout)
    filename = join(dirName, '%s.html' % readId)
    plotly.offline.plot(fig, show_link=False, filename=filename,
                        auto_open=interactive)


def plotEvaluation(y, x, names, yAxisLabel, title):
    """
    Plot the evaluation for different subsets of substrings.

    @param x: a C{list} of x coordinates.
    @param y: a C{list} of y coordinates.
    @param names: a C{list} of x tick names.
    @param yAxisLabel: a C{str} label for the y axis.
    @param title: a C{str} label for the title.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.grid()
    gridlines = ax.get_xgridlines() + ax.get_ygridlines()
    for line in gridlines:
        line.set_linestyle('-')
        line.set_color('lightgrey')
    ax.set_axisbelow(True)

    ax.plot(x, y, 'o',
            markerfacecolor='darkcyan', markeredgecolor='white',
            markersize=16)

    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(-0.5, len(names) - 0.5)
    ax.set_ylabel(yAxisLabel, fontsize=15)
    ax.set_title(title, fontsize=20)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    ax.set_xticks(range(8))
    ax.set_xticklabels(names, rotation=90, fontsize=15)
    plt.tick_params(axis='x', which='both', bottom='on', top='off',
                    labelbottom='on')
    plt.tick_params(axis='y', which='both', left='on', right='off',
                    labelbottom='on')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')


def plotEvaluations(fileList, title, totalTpr=True, fractionOfPdbCovered=True,
                    fractionOfStructuresCovered=True, correlation=True):
    """
    Plot the evaluation for different subsets.

    @param fileList: a C{list} of C{lists}, where each list contains the
        following: a C{str} filename with substrings that should be evaluated,
        a C{str} filename with substrings that should be evaluated as well as
        their rates of true and false positives, a C{str} filename of the pdb
        ss.txt file, a C{str} filename of the structure strings extracted from
        the pdbFile, a C{str} name of the structure type that should be
        considered. Must be one of 'AlphaHelix', 'AlphaHelix_3_10',
        'AlphaHelix_pi', 'ExtendedStrand' and a C{str} name of this subset of
        substrings.
    @param title: a C{str} title for the plot.
    @param totalTpr: if C{True}, evaluation about the total true positive rate
        will be performed.
    @param fractionOfPdbCovered: if C{True}, evaluation about the fraction of
        pdb structure sequences matched at least once will be performed.
    @param fractionOfStructuresCovered: if C{true}, evaluation about the
        fraction of known structures matched at least once will be performed.
    @param correlation: if C{True}, evaluation about the correlation between
        light matter scores of the known PDB finder and it's AC equivalent will
        be performed.
    """
    if totalTpr:
        x = []
        names = []
        for files in fileList:
            names.append(files[5])
            evaluation = PdbSubsetStatistics(files[0], files[1], files[2],
                                             files[3], files[4])
            x.append(evaluation.getTotalTpr())
        plotEvaluation(x, range(len(x)), names, 'Total TPR', title)

    if fractionOfPdbCovered:
        x = []
        names = []
        for files in fileList:
            names.append(files[5])
            evaluation = PdbSubsetStatistics(files[0], files[1], files[2],
                                             files[3], files[4])
            x.append(evaluation.getFractionOfPdbCovered())
        plotEvaluation(x, range(len(x)), names,
                       'Fraction of structures in PDB matched by a substring',
                       title)

    if fractionOfStructuresCovered:
        x = []
        names = []
        for files in fileList:
            names.append(files[5])
            evaluation = PdbSubsetStatistics(files[0], files[1], files[2],
                                             files[3], files[4])
            x.append(evaluation.getFractionOfStructuresCovered())
        plotEvaluation(x, range(len(x)), names,
                       ('Fraction of unique known %s matched by a substring' %
                        files[4]), title)

    if correlation:
        x = []
        names = []
        for files in fileList:
            names.append(files[5])
            evaluation = PdbSubsetStatistics(files[0], files[1], files[2],
                                             files[3], files[4])
            x.append(evaluation.getCorrelation())
        plotEvaluation(x, sum([4 * [i] for i in range(len(x))], []), names,
                       'Correlation coefficient between PDB LM scores and '
                       'scores using subset', title)
