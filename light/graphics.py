from __future__ import division

from warnings import warn
from random import uniform
import matplotlib.pyplot as plt
from matplotlib import patches, gridspec
import numpy as np
from operator import attrgetter
from textwrap import fill
from collections import defaultdict
from math import log10
from itertools import repeat

from light.backend import Backend
from light.colors import colors
from light.database import DatabaseSpecifier
from light.features import Landmark
from light.landmarks import ALL_LANDMARK_CLASSES
from light.parameters import FindParameters
from light.performance.overlap import CalculateOverlap
from light.performance import affinity
from light.bin_score import ALL_BIN_SCORE_CLASSES
from light.string import MultilineString
from light.significance import (
    Always, HashFraction, MaxBinHeight, MeanBinHeight)
from light.trig import ALL_TRIG_CLASSES

from dark.dimension import dimensionalIterator
from dark.fasta import FastaReads
from dark.reads import AAReadWithX

# Number of characters at which to wrap long lines in titles.
FILL_WIDTH = 120

ALL_FEATURES = [
    (feature.SYMBOL, feature.NAME) for feature in
    sorted(ALL_LANDMARK_CLASSES + ALL_TRIG_CLASSES,
           key=attrgetter('NAME'))]

# From http://www.randalolson.com/2014/06/28/how-to-make-beautiful-data-\
# visualizations-in-python-with-matplotlib/
#
# These are the "Tableau 20" colors as RGB.
TABLEAU20 = [
    (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the above RGB values to the [0, 1] range, the format matplotlib
# accepts.
for i in range(len(TABLEAU20)):
    r, g, b = TABLEAU20[i]
    TABLEAU20[i] = (r / 255.0, g / 255.0, b / 255.0)

# Make sure we have enough colors for our set of features.
assert len(ALL_FEATURES) <= len(TABLEAU20)

# The following maps from feature symbol to matplotlib color.
COLORS = dict(zip([feature[0] for feature in ALL_FEATURES], TABLEAU20))


def legendHandles(names):
    """
    Make a set of handles for a legend showing feature names and colors.

    @param names: A C{set} of C{light.features.Feature} instance names
        that should appear in the legend.
    @return: A C{list} of matplotlib C{Patch} instances for a legend.
    """
    return [patches.Patch(color=COLORS[symbol], label=name)
            for (symbol, name) in ALL_FEATURES if name in names]


def _rectangularPanel(rows, cols, title, makeSubPlot, equalizeXAxes=False,
                      equalizeYAxes=False, includeUpper=True,
                      includeLower=True, includeDiagonal=True, saveAs=False,
                      showFigure=True):
    """
    Create a rectangular panel of plots.

    @param rows: The C{int} number of rows that will appear in the panel.
    @param cols: The C{int} number of columns that will appear in the panel.
    @param title: A C{str} title for the figure.
    @param makeSubPlot: A function that accepts three arguments (an C{int} row,
        an C{int} column, and a matplotlib axes instance), which creates a
        sub-plot using the passed axes, and which returns a C{dict}
        containing information ('maxX', 'minX', 'maxY', 'minY', and an optional
        'title') about the created plot or C{None} if no plot should be shown
        in the given (row, col) location.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param equalizeYAxes: if C{True}, adjust the Y axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param includeUpper: If C{True} make calls to makeSubPlot for the
        upper triangle part of the C{rows} by C{cols} rectangle, when the row
        index is less than the column index. (Probably only useful for a square
        panel.) If C{False}, those sub-plots will be empty.
    @param includeLower: Like includeUpper, but for the lower triangle.
    @param includeDiagonal: Like includeUpper and includeLower, but for
        sub-plots on the diagonal.
    @param saveAs: A C{str} name for the file that the figure will be saved to.
    @param showFigure: If C{True} show the figure.
    """
    figure, ax = plt.subplots(rows, cols, squeeze=False, figsize=(150, 100))
    subplots = {}

    for row, col in dimensionalIterator((rows, cols)):
        if ((row < col and not includeUpper) or
                (row > col and not includeLower) or
                (row == col) and not includeDiagonal):
            subplots[(row, col)] = None
        else:
            subplots[(row, col)] = makeSubPlot(row, col, ax[row][col])

    if equalizeXAxes or equalizeYAxes:
        nonEmpty = [x for x in iter(subplots.values()) if x]
        title += '\n'
        if equalizeXAxes:
            maxX = max(subplot['maxX'] for subplot in nonEmpty)
            minX = min(subplot['minX'] for subplot in nonEmpty)
            title += 'X range: %s to %s' % (minX, maxX)
            if equalizeYAxes:
                title += ', '
        if equalizeYAxes:
            maxY = max(subplot['maxY'] for subplot in nonEmpty)
            minY = min(subplot['minY'] for subplot in nonEmpty)
            title += 'Y range: %s to %s' % (minY, maxY)

    # Post-process graphs to adjust axes, etc.
    for (row, col), subplot in subplots.items():
        a = ax[row][col]
        if subplot:
            try:
                subTitle = subplots[(row, col)]['title']
            except KeyError:
                # No title, no problem.
                pass
            else:
                a.set_title(fill(subTitle, 50), fontsize=10)
            if equalizeXAxes:
                a.set_xlim([minX, maxX])
                a.set_xticks([])
            if equalizeYAxes:
                a.set_ylim([minY, maxY])
                a.set_yticks([])
        else:
            # This subplot is not displayed.
            a.axis('off')

    figure.suptitle(title, fontsize=20)
    if saveAs:
        figure.set_size_inches(150, 100, forward=True, dpi=(300))
        figure.savefig(saveAs)
    if showFigure:
        figure.show()


def plotHistogramPanel(sequences, equalizeXAxes=True, equalizeYAxes=False,
                       findParams=None, showUpper=True, showLower=False,
                       showDiagonal=True, showMean=False, showMedian=False,
                       showStdev=False, showSignificanceCutoff=False,
                       showSignificantBins=True, saveAs=False, showFigure=True,
                       **kwargs):
    """
    Plot a square panel of histograms of matching hash offset deltas between
    all pairs of passed sequences.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param equalizeYAxes: if C{True}, adjust the Y axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param showUpper: If C{True}, show the sub-plots in the upper triangle.
        of the panel.
    @param showLower: If C{True}, show the sub-plots in the lower triangle.
        of the panel.
    @param showDiagonal: If C{True}, show the sub-plots on the diagonal.
        of the panel.
    @param showMean: If C{True} the mean will be plotted in red.
    @param showMedian: If C{True} the median will be plotted in orange.
    @param showStdev: If C{True} the standard deviation from the mean will be
        plotted in magenta.
    @param showSignificanceCutoff: If C{True} the significanceCutoff will be
        plotted in green.
    @param showSignificantBins: If C{True} the significant bins will be plotted
        in red.
    @param saveAs: A C{str} name for the file that the figure will be saved to.
    @param showFigure: If C{True} show the figure.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    """
    if isinstance(sequences, str):
        reads = list(FastaReads(sequences, readClass=AAReadWithX,
                     upperCase=True))
    else:
        reads = list(sequences)
    nReads = len(reads)
    # Make a new database. For now we don't allow an existing database to
    # be passed as we're going to use subject indices from 0 to nReads-1.
    # This shortcoming can be removed later.
    specifier = DatabaseSpecifier(allowInMemory=False)
    database = specifier.getDatabaseFromKeywords(subjects=reads, **kwargs)

    def makeSubPlot(row, col, ax):
        """
        @param row: The C{int} panel row.
        @param col: The C{int} panel column.
        @param ax: The matplotlib axis for the sub-plot.
        """
        return plotHistogram(reads[row], reads[col],
                             findParams=findParams, readsAx=ax,
                             showMean=showMean, showMedian=showMedian,
                             showStdev=showStdev,
                             showSignificanceCutoff=showSignificanceCutoff,
                             showSignificantBins=showSignificantBins,
                             database=database)

    return _rectangularPanel(
        nReads, nReads, 'Histogram panel', makeSubPlot,
        equalizeXAxes=equalizeXAxes, equalizeYAxes=equalizeYAxes,
        includeUpper=showUpper, includeLower=showLower,
        includeDiagonal=showDiagonal, saveAs=saveAs, showFigure=showFigure)


def plotHistogram(query, subject, findParams=None, readsAx=None,
                  showMean=False, showMedian=False, showStdev=False,
                  showSignificanceCutoff=False, showSignificantBins=False,
                  **kwargs):
    """
    Plot a histogram of matching hash offset deltas between a query and a
    subject.

    @param query: an AAReadWithX instance of the sequence of the query.
    @param subject: either an AAReadWithX instance of the sequence of the
        subject or a C{str} subject index in the database.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param showMean: If C{True} the mean will be plotted in red.
    @param showMedian: If C{True} the median will be plotted in orange.
    @param showStdev: If C{True} the standard deviation from the mean will be
        plotted in magenta.
    @param showSignificanceCutoff: If C{True} the significanceCutoff will be
        plotted in green.
    @param showSignificantBins: If C{True} the significant bins will be plotted
        in red.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    @raises IndexError: If C{subject} is an C{int} and that index does not
        exist in the database.
    """
    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)

    if isinstance(subject, str):
        subjectIndex = subject
        subject = database.getSubjectByIndex(subjectIndex)
    else:
        _, subjectIndex, subjectHashCount = database.addSubject(subject)

    result = database.find(query, findParams, storeFullAnalysis=True)

    try:
        histogram = result.analysis[subjectIndex]['histogram']
    except KeyError:
        print('Query %r and subject %r had no hashes in common.' % (
            query.id, subject.id))
    else:
        counts = [len(bin) for bin in histogram.bins]
        nBins = len(histogram.bins)
        centers = np.linspace(histogram.min, histogram.max, nBins,
                              endpoint=False)
        title = fill('%s vs %s' % (query.id, subject.id), FILL_WIDTH)

        if readsAx is None:
            fig = plt.figure()
            readsAx = fig.add_subplot(111)
            readsAx.set_title(title, fontsize=30)
            readsAx.set_xlabel('Offset delta (subject - query)', fontsize=14)
            readsAx.xaxis.tick_bottom()

        readsAx.vlines(centers, list(repeat(0, len(counts))), counts,
                       color='blue', linewidth=2)

        if showSignificantBins:
            significanceMethod = findParams.significanceMethod
            minHashCount = min(result.queryHashCount, subjectHashCount)

            if significanceMethod == 'Always':
                significance = Always()
            elif significanceMethod == 'HashFraction':
                significance = HashFraction(
                    histogram, minHashCount, findParams.significanceFraction)
            elif significanceMethod == 'MaxBinHeight':
                significance = MaxBinHeight(histogram, query, result.connector)
            elif significanceMethod == 'MeanBinHeight':
                significance = MeanBinHeight(histogram, query,
                                             result.connector)
            else:
                raise ValueError('Unknown significance method %r' %
                                 significanceMethod)

            for binIndex, bin_ in enumerate(histogram.bins):
                if significance.isSignificant(binIndex):
                    readsAx.vlines(centers[binIndex], 0, len(bin_),
                                   color='red', linewidth=2)

        mean = np.mean(counts)
        if showMean:
            readsAx.plot([histogram.min, histogram.max], [mean, mean], '-',
                         c='red')
        if showMedian:
            median = np.median(counts)
            readsAx.plot([histogram.min, histogram.max], [median, median], '-',
                         c='orange')
        if showStdev:
            stdev = np.std(counts)
            readsAx.plot([histogram.min, histogram.max],
                         [mean + stdev, mean + stdev], '-', c='magenta')

        if showSignificanceCutoff:
            sigAnalysis = result.analysis[subjectIndex]['significanceAnalysis']
            cutoff = sigAnalysis['significanceCutoff']
            readsAx.plot([histogram.min, histogram.max],
                         [cutoff, cutoff], '-', c='green')

        return {
            'minX': histogram.min,
            'maxX': histogram.max,
            'minY': 0,
            'maxY': max(counts),
            'result': result,
            'title': title,
        }


def plotHistogramLinePanel(sequences, equalizeXAxes=True, equalizeYAxes=False,
                           findParams=None, showUpper=True, showLower=False,
                           showDiagonal=True, showMean=False, showMedian=False,
                           showStdev=False, showSignificanceCutoff=False,
                           saveAs=False, showFigure=True, **kwargs):
    """
    Plot a square panel of histogram line plots of matching hash offset deltas
    between all pairs of passed sequences.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param equalizeYAxes: if C{True}, adjust the Y axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param showUpper: If C{True}, show the sub-plots in the upper triangle.
        of the panel.
    @param showLower: If C{True}, show the sub-plots in the lower triangle.
        of the panel.
    @param showDiagonal: If C{True}, show the sub-plots on the diagonal.
        of the panel.
    @param showMean: If C{True} the mean will be plotted in red.
    @param showMedian: If C{True} the median will be plotted in orange.
    @param showStdev: If C{True} the standard deviation from the mean will be
        plotted in magenta.
    @param showSignificanceCutoff: If C{True} the significanceCutoff will be
        plotted in green.
    @param saveAs: A C{str} name for the file that the figure will be saved to.
    @param showFigure: If C{True} show the figure.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    """
    if isinstance(sequences, str):
        reads = list(FastaReads(sequences, readClass=AAReadWithX,
                     upperCase=True))
    else:
        reads = list(sequences)
    nReads = len(reads)
    # Make a new database. For now we don't allow an existing database to
    # be passed as we're going to use subject indices from 0 to nReads-1.
    # This shortcoming can be removed later.
    specifier = DatabaseSpecifier(allowInMemory=False)
    database = specifier.getDatabaseFromKeywords(subjects=reads, **kwargs)

    def makeSubPlot(row, col, ax):
        """
        @param row: The C{int} panel row.
        @param col: The C{int} panel column.
        @param ax: The matplotlib axis for the sub-plot.
        """
        return plotHistogramLine(reads[row], reads[col],
                                 findParams=findParams,
                                 readsAx=ax, showMean=showMean,
                                 showMedian=showMedian, showStdev=showStdev,
                                 showSignificanceCutoff=showSignificanceCutoff,
                                 database=database)

    return _rectangularPanel(
        nReads, nReads, 'Histogram line panel', makeSubPlot,
        equalizeXAxes=equalizeXAxes, equalizeYAxes=equalizeYAxes,
        includeUpper=showUpper, includeLower=showLower,
        includeDiagonal=showDiagonal, saveAs=saveAs, showFigure=showFigure)


def plotHistogramLine(query, subject, findParams=None, readsAx=None,
                      showMean=False, showMedian=False, showStdev=False,
                      showSignificanceCutoff=False, **kwargs):
    """
    Plot a line where the height corresponds to the number of hashes in a
    histogram bin, but sorted by height.

    @param query: an AAReadWithX instance of the sequence of the query.
    @param subject: either an AAReadWithX instance of the sequence of the
        subject or an C{int} subject index in the database.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param showMean: If C{True} the mean will be plotted in red.
    @param showMedian: If C{True} the median will be plotted in orange.
    @param showStdev: If C{True} the standard deviation from the mean will be
        plotted in magenta.
    @param showSignificanceCutoff: If C{True} the significanceCutoff will be
        plotted in green.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    @raises IndexError: If C{subject} is an C{int} and that index does not
        exist in the database.
    """
    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)

    _, subjectIndex, _ = database.addSubject(subject)

    result = database.find(query, findParams, storeFullAnalysis=True)

    try:
        histogram = result.analysis[subjectIndex]['histogram']
    except KeyError:
        print('Query %r and subject %r had no hashes in common.' % (
            query.id, subject.id))
    else:
        counts = sorted(len(bin) for bin in histogram.bins)
        title = fill('%s vs %s' % (query.id, subject.id), FILL_WIDTH)

        if readsAx is None:
            fig = plt.figure()
            readsAx = fig.add_subplot(111)
            readsAx.set_title(title, fontsize=17)
            readsAx.set_ylabel('Number of hashes', fontsize=14)
            readsAx.xaxis.tick_bottom()

        readsAx.plot(range(len(counts)), counts)
        mean = np.mean(counts)
        if showMean:
            readsAx.plot([0, len(counts)], [mean, mean], '-',
                         c='red')
        if showMedian:
            median = np.median(counts)
            readsAx.plot([0, len(counts)], [median, median], '-',
                         c='orange')
        if showStdev:
            stdev = np.std(counts)
            readsAx.plot([0, len(counts)], [mean + stdev, mean + stdev],
                         '-', c='magenta')

        if showSignificanceCutoff:
            sigAnalysis = result.analysis[subjectIndex]['significanceAnalysis']
            cutoff = sigAnalysis['significanceCutoff']
            readsAx.plot([histogram.min, histogram.max],
                         [cutoff, cutoff], '-', c='green')

        return {
            'minX': 0,
            'maxX': len(counts),
            'minY': 0,
            'maxY': max(counts),
            'result': result,
            'title': title,
        }


def plotHistogramLines(sequences, findParams=None, **kwargs):
    """
    Plot lines where the height corresponds to the number of hashes in a
    histogram bin, but sorted by height, do that for many reads.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    fig = plt.figure(figsize=(15, 10))
    readsAx = fig.add_subplot(111)
    if isinstance(sequences, str):
        reads = list(FastaReads(sequences, readClass=AAReadWithX,
                     upperCase=True))
    else:
        reads = list(sequences)
    # Make a new database. For now we don't allow an existing database to
    # be passed as we're going to use subject indices from 0 to nReads-1.
    # This shortcoming can be removed later.
    specifier = DatabaseSpecifier(allowInMemory=False)
    database = specifier.getDatabaseFromKeywords(subjects=reads, **kwargs)

    for read in reads:
        result = database.find(read, findParams, storeFullAnalysis=True)
        for subjectIndex in map(str, range(len(reads))):
            try:
                analysis = result.analysis[subjectIndex]
            except KeyError:
                subject = database.getSubjectByIndex(subjectIndex)
                print('Query %r and subject %r had no hashes in common.' % (
                    read.id, subject.id))
            else:
                histogram = analysis['histogram']
                counts = sorted([len(bin) for bin in histogram.bins])
                readsAx.plot(range(len(counts)), counts)

    readsAx.set_title('Histogram line plot', fontsize=17)
    readsAx.set_ylabel('Number of hashes', fontsize=14)
    readsAx.xaxis.tick_bottom()


def plotFeatureSquare(read, findParams=None, readsAx=None, **kwargs):
    """
    Plot the positions of landmark and trigpoint pairs on a sequence in a
    square.

    @param read: A C{dark.reads.Read} instance.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    fig = plt.figure(figsize=(15, 15))
    readsAx = readsAx or fig.add_subplot(111)

    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    backend = Backend()
    backend.configure(database.dbParams)
    result = database.find(read, findParams, storeFullAnalysis=True)
    scannedQuery = backend.scan(result.query)

    # Plot a light grey diagonal line, bottom left to top right.
    readsAx.plot([0, len(read.sequence)], [0, len(read.sequence)], '-',
                 color='#f0f0f0', linewidth=1)

    scatterX = []
    scatterY = []
    scatterColors = []
    namesSeen = set()
    landmarks = set()

    for landmark, trigPoint in backend.getScannedPairs(scannedQuery):
        # Add jitter to the Y offset so we can see more trig points that
        # occur at the same offset.
        scatterX.append(landmark.offset)
        scatterY.append(trigPoint.offset + uniform(-0.4, 0.4))
        if isinstance(trigPoint, Landmark):
            # Color the trig points that are actually landmarks a light grey.
            scatterColors.append((0.7, 0.7, 0.7))
        else:
            scatterColors.append(COLORS[trigPoint.symbol])
        landmarks.add((landmark.symbol, landmark.offset, landmark.length))
        namesSeen.update([landmark.name, trigPoint.name])

    # Display all landmarks.
    for symbol, offset, length in landmarks:
        # Add jitter to the landmark offset so we can see landmarks that
        # might otherwise overlap.
        offset += uniform(-0.4, 0.4)
        readsAx.plot([offset, offset + length - 1], [offset, offset], '-',
                     color=COLORS[symbol], linewidth=2)

    # Show trig points for all landmarks using a scatter plot.
    readsAx.scatter(scatterX, scatterY, c=scatterColors, edgecolors='none')

    # Set labels, titles, axis limits, legend, etc.
    totalCoveredResidues = len(scannedQuery.coveredIndices())
    readsAx.set_title(
        '%s\n Length: %d, covered residues: %s' % (
            fill(read.id, FILL_WIDTH), len(read), totalCoveredResidues),
        fontsize=17)
    readsAx.set_xlabel('Offset (landmarks)', fontsize=14)
    readsAx.set_ylabel('Offset (trig points)', fontsize=14)
    readsAx.set_xlim(-0.2, len(read.sequence) + 0.2)
    readsAx.set_ylim(0, len(read.sequence))
    if namesSeen:
        readsAx.legend(handles=legendHandles(namesSeen), loc=2)
    readsAx.grid()

    return scannedQuery


def plotHorizontalPairPanel(sequences, findParams=None, equalizeXAxes=True,
                            showSignificant=True, showInsignificant=True,
                            showBestBinOnly=False, showUpper=True,
                            showLower=False, showDiagonal=True, saveAs=False,
                            showFigure=True, **kwargs):
    """
    Plot a panel of paired horizontally aligned sequences showing matching
    hashes. The individual sub-plots are produced by
    PlotHashesInSubjectAndRead.plotHorizontal.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param findParams: An instance of C{FindParameters} or C{None}, to use the
        default find parameters.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param showSignificant: If C{True}, hashes from significant bins will
        be included in the set of hashes that match query and subject.
    @param showInsignificant: If C{True}, hashes from insignificant bins will
        be included in the set of hashes that match query and subject.
    @param showBestBinOnly: If C{True}, only show the bin with the best
        score. Warn if there are multiple bins with the same high score.
    @param showUpper: If C{True}, show the sub-plots in the upper triangle.
        of the panel.
    @param showLower: If C{True}, show the sub-plots in the lower triangle.
        of the panel.
    @param showDiagonal: If C{True}, show the sub-plots on the diagonal.
        of the panel.
    @param saveAs: A C{str} name for the file that the figure will be saved to.
    @param showFigure: If C{True} show the figure.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    """
    if isinstance(sequences, str):
        reads = list(FastaReads(sequences, readClass=AAReadWithX,
                     upperCase=True))
    else:
        reads = list(sequences)
    nReads = len(reads)
    # Make a new database. For now we don't allow an existing database to
    # be passed as we're going to use subject indices from 0 to nReads-1.
    # This shortcoming can be removed later.
    specifier = DatabaseSpecifier(allowInMemory=False)
    database = specifier.getDatabaseFromKeywords(subjects=reads, **kwargs)

    def makeSubPlot(row, col, ax):
        """
        @param row: The C{int} panel row.
        @param col: The C{int} panel column.
        @param ax: The matplotlib axis for the sub-plot.
        """
        plotter = PlotHashesInSubjectAndRead(
            reads[row], reads[col], findParams,
            showSignificant=showSignificant,
            showInsignificant=showInsignificant,
            showBestBinOnly=showBestBinOnly, database=database)
        return plotter.plotHorizontal(ax)

    return _rectangularPanel(
        nReads, nReads, 'Horizontal pairs panel', makeSubPlot,
        equalizeXAxes=equalizeXAxes, equalizeYAxes=False,
        includeUpper=showUpper, includeLower=showLower,
        includeDiagonal=showDiagonal, saveAs=saveAs, showFigure=showFigure)


class PlotHashesInSubjectAndRead(object):
    """
    A class which plots visualisations of the hashes in subject and query.
    It collects three types of hashes: 1) Hashes in the query that don't
    match in the subject. 2) Hashes in the subject that don't match in the
    query. 3) Hashes that match in subject and query. These can subsequently
    be plotted.

    @param query: An AAReadWithX instance of the sequence of the query.
    @param subject: An AAReadWithX instance of the sequence of the subject.
    @param findParams: An instance of C{FindParameters} or C{None}, to use the
        default find parameters.
    @param showSignificant: If C{True}, hashes from significant bins will
        be included in the set of hashes that match query and subject.
    @param showInsignificant: If C{True}, hashes from insignificant bins will
        be included in the set of hashes that match query and subject.
    @param showBestBinOnly: If C{True}, only show the bin with the best
        score. Warn if there are multiple bins with the same high score.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    def __init__(self, query, subject, findParams=None, showSignificant=True,
                 showInsignificant=True, showBestBinOnly=False, **kwargs):
        self.query = query
        self.subject = subject

        database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        _, subjectIndex, _ = database.addSubject(subject)
        findParams = findParams or FindParameters()
        self.result = database.find(self.query, findParams,
                                    storeFullAnalysis=True)

        # Use an 'in' test here, as self.result.analysis is a defaultdict.
        # Relying on a KeyError via directly accessing
        # self.result.analysis[subjectIndex] won't work as the key would
        # just be created and no KeyError occurs.
        if subjectIndex in self.result.analysis:
            analysis = self.result.analysis[subjectIndex]
            self.score = analysis['bestBinScore']
            self.significantBinCount = len(analysis['significantBins'])

            # If showBestBinOnly is true, we need significantBins to be
            # non-empty in order to have a bin to examine.  This is because
            # insignificant bins do not have a score computed for them.
            # Note that the scores of the bins are sorted (best first) in
            # the Result class, so the first bin is the one with the best
            # score.
            if showBestBinOnly and analysis['significantBins']:
                self.bins = [analysis['significantBins'][0]['bin']]
                if len(analysis['significantBins']) > 1 and (
                        analysis['significantBins'][0]['score'] ==
                        analysis['significantBins'][1]['score']):
                    warn('Multiple bins share the best score (%f). '
                         'Displaying just one of them.' %
                         analysis['significantBins'][0]['score'],
                         RuntimeWarning)
            elif showSignificant:
                if showInsignificant:
                    self.bins = [bin_ for bin_ in analysis['histogram'].bins]
                else:
                    self.bins = [bin_['bin'] for bin_ in
                                 analysis['significantBins']]
            elif showInsignificant and not showSignificant:
                significantBinIndices = set(
                    [bin_['index'] for bin_ in analysis['significantBins']])
                self.bins = [bin_ for i, bin_ in
                             enumerate(analysis['histogram'].bins) if i not
                             in significantBinIndices]

        else:
            # The subject was not matched.
            self.score = 0.0
            self.significantBinCount = 0
            self.bins = []

        self.queryHashes = self.result.nonMatchingHashes
        backend = Backend()
        backend.configure(database.dbParams)
        self.subjectHashes = backend.getHashes(backend.scan(subject))
        self.matchingHashes = defaultdict(list)

        for bin_ in self.bins:
            for match in bin_:
                hash_ = backend.hash(match['subjectLandmark'],
                                     match['subjectTrigPoint'])
                self.matchingHashes[hash_].append(match)
                try:
                    del self.subjectHashes[hash_]
                except KeyError:
                    continue

    def plotGraph(self, readsAx=None):
        """
        Plots the hashes in a subject and a query. Matching hashes are plotted
        as a scatter plot, with hashes from each histogram bin coloured
        separately. Non-matching hashes in the subject are plotted on the
        x-axis, non matching hashes in the query are plotted on the y-axis.

        @param readsAx: If not C{None}, use this as the subplot for displaying
            reads.
        """
        height = (len(self.query) * 15) / len(self.subject)
        fig = plt.figure(figsize=(15, height))
        readsAx = readsAx or fig.add_subplot(111)

        for hashInfoList in self.queryHashes.values():
            for hashInfo in hashInfoList:
                landmarkOffset = hashInfo[0].offset
                readsAx.plot(landmarkOffset + uniform(-0.4, 0.4), 0, 'o',
                             markerfacecolor='black', markeredgecolor='white')

        for hashInfoList in self.subjectHashes.values():
            for hashInfo in hashInfoList:
                landmarkOffset = hashInfo[0].offset
                readsAx.plot(0, landmarkOffset + uniform(-0.4, 0.4), 'o',
                             markerfacecolor='black', markeredgecolor='white')

        nonEmptyBins = [b for b in self.bins if len(b)]
        binColors = colors.color_palette('hls', len(nonEmptyBins))
        for bin_, binColor in zip(nonEmptyBins, binColors):
            for match in bin_:
                readsAx.plot(
                    match['subjectLandmark'].offset + uniform(-0.4, 0.4),
                    match['queryLandmark'].offset + uniform(-0.4, 0.4),
                    'o', markerfacecolor=binColor, markeredgecolor='white')

        firstTitleLine = fill('Hashes from matching %s against %s' % (
            self.query.id, self.subject.id), FILL_WIDTH)
        readsAx.set_title(
            '%s\nScore: %.4f, HSPs: %d\n'
            'Hashes: matching=%d, subject-only=%d, query-only=%d' % (
                firstTitleLine, self.score or 0.0, self.significantBinCount,
                len(self.matchingHashes), len(self.subjectHashes),
                len(self.queryHashes)), fontsize=17)
        readsAx.set_ylabel('Query: %s' % self.query.id, fontsize=14)
        readsAx.set_xlabel('Subject: %s' % self.subject.id, fontsize=14)
        readsAx.set_xlim(-0.5, len(self.subject))
        readsAx.set_ylim(-0.5, len(self.query))
        readsAx.grid()

    def plotHorizontal(self, ax=None, addJitter=True):
        """
        Plot a graph with the subject and query drawn as two horizontal lines.

        @param ax: If not C{None}, use this as the subplot.
        @param addJitter: If C{True}, add a small amount of jitter to the
            x-axis offset used to plot landmarks and trig points to decrease
            the probability that plotted lines completely overlap one another.
        """

        # In the code below, the following abbreviations are used in
        # variable names:
        #
        #   qy: query
        #   sj: subject
        #   lm: landmark
        #   tp: trig point
        #   diag: diagonal
        #
        # For consistency, query handling is always done before subject
        # handling, and landmarks are handled before trig points.

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            createdAx = True
        else:
            createdAx = False
        namesSeen = set()
        qyY = 0.0
        sjY = 1.0
        tpY = 0.015
        missY = 0.1
        maxLen = max([len(self.query), len(self.subject)])
        verticalPad = 0.1
        horizontalPad = 0.02 * maxLen

        # A line for query hashes that match the subject.
        ax.plot([0.0, len(self.query)], [qyY, qyY], '-', color='#bbbbbb')

        # A lower line for query hashes that don't match the subject.
        ax.plot([0.0, len(self.query)], [qyY - missY, qyY - missY],
                '-', color='#cccccc')

        # A line for subject hashes that match the query.
        ax.plot([0.0, len(self.subject)], [sjY, sjY], '-', color='#bbbbbb')

        # A higher line for subject hashes that don't match the query.
        ax.plot([0.0, len(self.subject)], [sjY + missY, sjY + missY],
                '-', color='#cccccc')

        # Keep track of features plotted in the query and subject, so we
        # don't plot them in the missed feature areas. Similar for the
        # diagonal lines showing matches.
        qyPlotted = set()
        sjPlotted = set()
        diagPlotted = set()

        def plotFeature(feature, y, queryOrSubject):
            """
            Plot a feature.

            @param feature: A C{light.features._Feature} instance.
            @param y: The C{float} Y offset.
            @param queryOrSubject: A C{str}, either 'query' or 'subject'.
            @raise KeyError: If C{queryOrSubject} is not 'query' or 'subject'.
            """
            seen = {
                'query': qyPlotted,
                'subject': sjPlotted,
            }[queryOrSubject]

            if feature in seen:
                return

            seen.add(feature)
            color = COLORS[feature.symbol]

            if isinstance(feature, Landmark):
                ax.plot([feature.offset, feature.offset + feature.length - 1],
                        [y, y], '-', color=color, linewidth=3)
            else:
                ax.plot([feature.offset, feature.offset], [y + tpY, y - tpY],
                        '-', color=color, linewidth=3)

        # Plot matching hashes on the query and subject horizontal lines,
        # and plot diagonal lines between the two (using colors corresponding
        # to histogram bins).
        #
        # Note that the dictionary keys of the items in histogram bins are
        # different from those in the query-only hashes and subject-only
        # hashes (processed later) as those are the result of calling
        # getHashes in a backend.

        nonEmptyBins = [b for b in self.bins if len(b)]
        binColors = colors.color_palette('hls', len(nonEmptyBins))
        jitterRange = 0.005 * max(len(self.query), len(self.subject))

        for bin_, binColor in zip(nonEmptyBins, binColors):
            for match in bin_:
                qlm = match['queryLandmark']
                qtp = match['queryTrigPoint']
                slm = match['subjectLandmark']
                stp = match['subjectTrigPoint']

                namesSeen.update([qlm.name, qtp.name])

                # Plot the landmark and trig point in the query and subject.
                plotFeature(qlm, qyY, 'query')
                plotFeature(qtp, qyY, 'query')
                plotFeature(slm, sjY, 'subject')
                plotFeature(stp, sjY, 'subject')

                # Add x-axis offset jitter to diagonal line plotting, so we
                # can see more of them in case of overlap.
                xJitter = (
                    uniform(-jitterRange, jitterRange) if addJitter else 0.0)

                # Dashed diagonal line connecting trig point in query and
                # subject.  Plot this before the diagonal landmark lines as
                # it's arguably less important (there's a chance we will
                # paint over it with the landmark lines).
                key = (qtp, stp, 'trigPoint')
                if key not in diagPlotted:
                    diagPlotted.add(key)
                    ax.plot([qtp.offset + xJitter, stp.offset + xJitter],
                            [qyY, sjY], '--', color=binColor)

                # Solid diagonal line connecting start of landmark in query
                # and subject.
                key = (qlm, slm, 'landmark')
                if key not in diagPlotted:
                    diagPlotted.add(key)
                    ax.plot([qlm.offset + xJitter, slm.offset + xJitter],
                            [qyY, sjY], '-', color=binColor)

        # Query-only hashes, plotted just below (-missY) the query line.
        #
        # Note that the keys of the items in queryHashes are different from
        # those in the histogram bins (processed above) as these hashes are
        # the result of calling getHashes in a backend.
        for hashInfoList in self.queryHashes.values():
            for hashInfo in hashInfoList:
                lm, tp = hashInfo
                namesSeen.update([lm.name, tp.name])
                plotFeature(lm, qyY - missY, 'query')
                plotFeature(tp, qyY - missY, 'query')

        # Subject-only hashes, plotted just above (+missY) the subject line.
        #
        # Note that, as with query-only hashes, the keys of the items in
        # subjectHashes are different from those in the histogram bins
        # (processed above) as these hashes are the result of calling
        # getHashes in a backend.
        for hashInfoList in self.subjectHashes.values():
            for hashInfo in hashInfoList:
                lm, tp = hashInfo
                namesSeen.update([lm.name, tp.name])
                plotFeature(lm, sjY + missY, 'subject')
                plotFeature(tp, sjY + missY, 'subject')

        if createdAx:
            if namesSeen:
                ax.legend(handles=legendHandles(namesSeen),
                          bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3,
                          ncol=2, borderaxespad=0.5)
            ax.set_xlabel(fill('%s (top) vs %s (bottom)' %
                               (self.subject.id, self.query.id)), fontsize=14)
        minX = -horizontalPad
        maxX = maxLen + horizontalPad
        ax.set_xlim(minX, maxX)

        minY = qyY - missY - verticalPad
        maxY = sjY + missY + verticalPad
        ax.set_ylim(minY, maxY)
        ax.set_yticks([])

        return {
            'minX': minX,
            'maxX': maxX,
            'minY': minY,
            'maxY': maxY,
            'title': fill('%s vs %s' % (self.subject.id[:20],
                                        self.query.id[:20])),
            'colors': binColors,
            'namesSeen': namesSeen,
        }


def plotLandmarksInSequences(sequences, maxTickLabelLength=None, **kwargs):
    """
    Plot the positions of landmarks and trig points on many sequences, with
    sequences stacked above each other.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param maxTickLabelLength: An C{int} limit on the length of tick labels
        on the Y axis. If C{None}, labels will not be truncated.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    if isinstance(sequences, str):
        reads = list(FastaReads(sequences, readClass=AAReadWithX,
                     upperCase=True))
    else:
        reads = list(sequences)

    nReads = len(reads)
    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    backend = Backend()
    backend.configure(db.dbParams)
    fig = plt.figure(figsize=(15, nReads / 3.0))
    ax = fig.add_subplot(111)
    namesSeen = set()
    maxLen = 0
    yticks = []

    for i, read in enumerate(reads):
        y = nReads - 1 - i
        yticks.append(read.id if maxTickLabelLength is None else
                      read.id[:maxTickLabelLength])
        readLen = len(read)
        if readLen > maxLen:
            maxLen = readLen
        plt.plot([0, len(read.sequence)], [y, y], '-', linewidth=0.5,
                 color='grey')
        scannedRead = backend.scan(read)
        # Landmarks are drawn as colored horizontal lines.
        for landmark in scannedRead.landmarks:
            namesSeen.add(landmark.name)
            plt.plot([landmark.offset, landmark.offset + landmark.length - 1],
                     [y, y], '-', color=COLORS[landmark.symbol], linewidth=2)
        # Trig points are drawn as small colored vertical lines.
        for trigPoint in scannedRead.trigPoints:
            namesSeen.add(trigPoint.name)
            plt.plot([trigPoint.offset, trigPoint.offset],
                     [y - 0.125, y + 0.125], '-',
                     color=COLORS[trigPoint.symbol], linewidth=2)

    ax.set_title('Landmarks and features\n', fontsize=17)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.xaxis.grid()
    ax.set_ylim(-0.1, nReads - 1 + 0.1)
    ax.set_xlim(0, maxLen)
    yticks.reverse()
    ax.set_yticklabels(yticks)
    ax.set_yticks(range(nReads))
    # Add a legend above left on the plot.
    if namesSeen:
        ax.legend(handles=legendHandles(namesSeen),
                  bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3, ncol=2,
                  borderaxespad=0.5)

    plt.tick_params(axis='x', which='both', bottom='off', top='off',
                    labelbottom='on')


def confusionMatrix(confusionMatrix):
    """
    Plot a confusion matrix (http://en.wikipedia.org/wiki/Confusion_matrix)
    to visualize the accuracy of the clustering.

    @param confusionMatrix: A confusion matrix, as returned from
        C{light.performance.cluster.ClusterAnalysis}.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    plt.imshow(confusionMatrix, interpolation='nearest', cmap=plt.cm.GnBu)
    plt.title('Confusion matrix \n')
    colorbar = plt.colorbar(ticks=[list(range(0,
                            max(confusionMatrix.max(axis=1)) + 1))])
    colorbar.outline.remove()
    plt.tick_params(axis='x', which='both', bottom='off', top='off',
                    labelbottom='off', labeltop='on')
    plt.tick_params(axis='y', which='both', left='off', right='off',
                    labelbottom='on')
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)


def featureComparison(ssAARead, print_=True, **kwargs):
    """
    A function which provides plotting options for sequences, given the
    predicted secondary structures from PDB and our features.

    Abbreviations in the pdb secondary structures:
    H = alpha-helix
    B = residue in isolated beta-bridge
    E = extended strand, participates in beta ladder
    G = 3-helix (310 helix)
    I = 5-helix (pi-helix)
    T = hydrogen bonded turn
    S = bend

    @param ssAARead: A C{light.performance.overlap.SSAARead} instance.
    @param print_: A C{bool} indicating if the result of the overlap
        calculation should be printed.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    ssStructureColors = {
        'H': (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
        'B': '#848484',
        'E': (0.596078431372549, 0.8745098039215686, 0.5411764705882353),
        'G': (0.6823529411764706, 0.7803921568627451, 0.9098039215686274),
        'I': (1.0, 0.4980392156862745, 0.054901960784313725),
        'T': '#848484',
        'S': '#848484',
    }

    ssStructure = ssAARead.structure

    # parse the ssStructure so that it can be plotted easily.
    all_ = defaultdict(list)

    previous = None
    start = 0
    for i, item in enumerate(ssStructure):
        # item is the last item of the sequence
        try:
            ssStructure[i + 1]
        except IndexError:
            if item == previous:
                all_[item].append([start, i])
            else:
                all_[item].append([i, i])
            break
        if item == previous and ssStructure[i + 1] != item:
            all_[item].append([start, i])
        # item is the only item of the sequence
        elif item != previous and ssStructure[i + 1] != item:
            previous = item
            all_[item].append([i, i])
        # item is the first one in the sequence:
        elif item != previous and ssStructure[i + 1] == item:
            previous = item
            start = i

    # set up database and scan read.
    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    backend = Backend()
    backend.configure(db.dbParams)
    scannedRead = backend.scan(ssAARead)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    lmNames = db.dbParams.landmarkFinderNames()
    tpNames = db.dbParams.trigPointFinderNames()
    yticks = ['S', 'T', 'E', 'B', 'I', 'G', 'H', ' '] + lmNames + tpNames
    ytickLabels = (['Bend', 'H-bonded turn', 'BetaStrand (?)', 'BetaBridge',
                    'AlphaHelixPi', 'AlphaHelix_3_10', 'AlphaHelix', ' '] +
                   lmNames + tpNames)
    title = ssAARead.id

    for i, item in enumerate(ytickLabels):
        plt.plot([0, len(ssAARead.sequence)], [i, i], '-', linewidth=0.5,
                 color='grey')

    for landmark in scannedRead.landmarks:
        y = yticks.index(landmark.name)
        if landmark.name != 'AminoAcidsLm':
            plt.plot([landmark.offset, landmark.offset + landmark.length - 1],
                     [y, y], '-', color=COLORS[landmark.symbol], linewidth=4)
        else:
            plt.plot([landmark.offset, landmark.offset + landmark.length],
                     [y - 0.125, y + 0.125], '-',
                     color=COLORS[landmark.symbol], linewidth=2)
    for trigPoint in scannedRead.trigPoints:
        y = yticks.index(trigPoint.name)
        plt.plot([trigPoint.offset, trigPoint.offset], [y - 0.125, y + 0.125],
                 '-', color=COLORS[trigPoint.symbol], linewidth=2)

    plt.plot([0, len(ssAARead.sequence)], [7, 7], '-', linewidth=3,
             color='black')

    for feature, offsets in all_.items():
        y = yticks.index(feature)
        if feature != ' ':
            for offset in offsets:
                plt.plot([offset[0], offset[1]], [y, y],
                         color=ssStructureColors[feature], linewidth=4)

    ax.set_yticklabels(ytickLabels)
    ax.set_yticks(range(len(yticks)))
    plt.title(title, fontsize=15)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.xaxis.grid()
    ax.set_ylim(-0.1, len(yticks) - 1 + 0.1)
    ax.set_xlim(0, len(ssAARead.sequence))

    if print_:
        aaSeqF, ssSeqF, intersects = CalculateOverlap.getFeatures(ssAARead,
                                                                  **kwargs)
        CalculateOverlap.calculateFraction(aaSeqF, ssSeqF, intersects,
                                           print_=True)


class SequenceFeatureAnalysis:
    """
    Analyze a sequence for its features and provide methods to plot and
    print useful summaries.

    @param sequence: Either a C{str} AA sequence or an C{AAReadWithX} instance.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    def __init__(self, sequence, **kwargs):
        if isinstance(sequence, str):
            sequence = AAReadWithX('id', sequence)

        self.sequence = sequence

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(db.dbParams)
        scannedSequence = backend.scan(sequence)

        self.offsets = defaultdict(set)
        self.landmarkNames = set()
        self.trigPointNames = set()

        for feature in scannedSequence.landmarks:
            self.landmarkNames.add(feature.name)
            self.offsets[feature.name].update(feature.coveredOffsets())

        for feature in scannedSequence.trigPoints:
            self.trigPointNames.add(feature.name)
            self.offsets[feature.name].update(feature.coveredOffsets())

        self.landmarksNotFound = (
            set(kwargs.get('landmarks', [])) - self.landmarkNames)

        self.trigPointsNotFound = (
            set(kwargs.get('trigPoints', [])) - self.trigPointNames)

        self.allNames = (sorted(self.landmarkNames) +
                         sorted(self.trigPointNames))

    def pairwiseOffsetOverlaps(self):
        """
        For each pair of features, work out the fraction of AAs that are
        covered by both features compared to the number covered by either.

        @return: A numpy array containing the overlap fractions. The rows
            and columns are in the order of names in self.allNames. The
            matrix will be symmetric, with 1.0 diagonal values.
        """
        allNames = self.allNames
        overlaps = np.zeros((len(allNames), len(allNames)))
        offsets = self.offsets

        for i, name1 in enumerate(allNames):
            for j, name2 in enumerate(allNames):
                if i == j:
                    overlaps[i, i] = 1.0
                else:
                    inCommon = len(offsets[name1] & offsets[name2])
                    inTotal = len(offsets[name1] | offsets[name2])
                    try:
                        overlaps[i, j] = inCommon / inTotal
                    except ZeroDivisionError:
                        pass

        return overlaps

    def plotPairwiseHeatmap(self):
        """
        Plot a heatmap showing the pairwise feature offset overlap.
        """
        overlaps = self.pairwiseOffsetOverlaps()
        # Set the diagonal to 0.0 to not wildly distort the color bar with
        # 1.0 values.
        for i in range(len(self.allNames)):
            overlaps[i][i] = 0.0

        img = plt.imshow(overlaps, interpolation='nearest', cmap=plt.cm.GnBu)

        left, right, bottom, top = img.get_extent()

        # Draw feature names on the left and bottom of the image.
        for i, name in enumerate(self.allNames):
            plt.text(left - 6, top + i + 0.75, name, fontsize=15)
            plt.text(left + i + 0.3, bottom + 1, name, fontsize=15,
                     rotation='vertical')

        plt.xticks([])
        plt.yticks([])

        # Draw lines to mark off landmarks from trig points.
        plt.axhline(y=bottom - len(self.trigPointNames), linestyle=':',
                    color='#888888')
        plt.axvline(x=right - len(self.trigPointNames), linestyle=':',
                    color='#888888')

        plt.colorbar(img)
        plt.title('Pairwise feature offset overlap\n', fontsize=18)
        plt.show()

    def printDensities(self, margin=''):
        """
        Produce a table of offset densities for each feature.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} summarizing feature densities and the number of
            offsets that are not covered by any feature.
        """
        result = MultilineString(margin=margin, indent='  ')
        result.append('Feature densities:')
        result.indent()

        totalOffsetCount = len(self.sequence)
        uncoveredOffsets = set(range(totalOffsetCount))
        unique = self.uniqueOffsets()
        maxName = max(len(name) for name in self.offsets)
        maxCount = int(log10(totalOffsetCount)) + 1
        maxOffsetCount = int(log10(max(
            len(offsets) for offsets in self.offsets.values()))) + 1
        maxUniqueOffsetCount = int(log10(max(
            len(offsets) for offsets in unique.values()))) + 1

        result.append('%*s OVERALL%*s   UNIQUE' % (
            maxName + 1, '',
            10 + maxOffsetCount + maxCount - len('OVERALL'), ''))
        for name in self.allNames:
            offsets = self.offsets[name]
            offsetCount = len(offsets)
            result.append('%-*s %5.2f%% (%*d/%d)   %5.2f%% (%*d/%d)' % (
                maxName + 1, name + ':',
                offsetCount / totalOffsetCount * 100.0,
                maxOffsetCount, offsetCount, totalOffsetCount,
                len(unique[name]) / totalOffsetCount * 100.0,
                maxUniqueOffsetCount, len(unique[name]), totalOffsetCount))
            uncoveredOffsets -= offsets

        uncoveredOffsetCount = len(uncoveredOffsets)

        result.outdent()
        result.append('')
        result.append(
            '%5.2f%% (%d/%d) of sequence offsets were not covered by any '
            'feature.' % (uncoveredOffsetCount / totalOffsetCount * 100.0,
                          uncoveredOffsetCount, totalOffsetCount))

        return str(result)

    def uniqueOffsets(self):
        """
        Calculate the set of offsets that each feature is unique in covering.

        @return: A C{dict} whose keys are feature names and whose values are
            sets of offsets that only that feature covers. There is a key
            present even for finders that found no features.
        """
        unique = {}
        offsets = self.offsets

        # First, copy the offsets for each feature.
        for name in self.allNames:
            unique[name] = set(offsets[name])

        # Then, from each finder, subtract the offsets of all other finders.
        for name in self.allNames:
            offsetsForThisFinder = unique[name]
            for otherName in self.allNames:
                if otherName != name:
                    offsetsForThisFinder -= offsets[otherName]

        # Add in an empty set for those finders that found nothing.
        for name in self.landmarksNotFound | self.trigPointsNotFound:
            unique[name] = set()

        return unique


def compareScores(subject, query, binScoreMethods=None,
                  plotHashesInSubjectAndRead=True, showBestBinOnly=True,
                  showHistogram=True, findParams=None, showInsignificant=False,
                  showFeatures=False, showScoreAnalysis=True,
                  showSignificantBinsDetails=False,
                  showBestBinFeatureInfo=False, **kwargs):
    """
    Plot the features in two sequences, the hashes in their match, and
    the scores for the match calculated under different score methods.

    @param subject: An C{AAReadWithX} instance.
    @param query: An C{AAReadWithX} instance.
    @param binScoreMethods: A C{list} of bin score methods to use, or C{None}
        to use all score methods.
    @param plotHashesInSubjectAndRead: If C{True}, plot the hashes in the
        subject and read showing the match.
    @param showBestBinOnly: If C{True} (and C{PlotHashesInSubjectAndRead} is
        also C{True}), only show details of the best bin in the match between
        query and subject.
    @param showHistogram: If C{True}, plot the delta offset histogram.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param showInsignificant: If C{True}, hashes from insignificant bins will
        be included in the set of hashes that match query and subject.
    @param showFeatures: If C{True}, show a separate plot of features in the
        sequences.
    @param showScoreAnalysis: If C{True}, print details of the score analysis.
    @param showSignificantBinsDetails: If C{True}, print information about the
        contents of the significant bins in the histogram.
    @param showBestBinFeatureInfo: If C{True} print which features are found
        in the best bin.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    _, subjectIndex, _ = db.addSubject(subject)
    binScoreMethods = binScoreMethods or sorted(
        cls.__name__ for cls in ALL_BIN_SCORE_CLASSES)
    findParams = findParams or FindParameters()

    # Calculate a score for each score method.
    for binScoreMethod in binScoreMethods:
        print('%s:' % binScoreMethod)
        findParams.binScoreMethod = binScoreMethod
        result = db.find(query, findParams, storeFullAnalysis=True).analysis
        if subjectIndex in result:
            significantBins = result[subjectIndex]['significantBins']
            if significantBins:
                if showScoreAnalysis:
                    analysis = significantBins[0]['scoreAnalysis']
                    print(analysis['scoreClass'].printAnalysis(analysis))
                else:
                    print('Score:', result[subjectIndex]['bestScore'])

                if showSignificantBinsDetails:
                    print('There are %d significant bins:\n  %s' % (
                        len(significantBins),
                        '\n  '.join(['index=%d score=%.4f binCount=%d' %
                                    ((b['index'], b['score'], len(b['bin'])))
                                    for b in significantBins])))
                if showBestBinFeatureInfo:
                    bestBinInfo = significantBins[0]['bin']
                    features = defaultdict(int)
                    for hash_ in bestBinInfo:
                        features[hash_['subjectLandmark'].name] += 1
                        features[hash_['subjectTrigPoint'].name] += 1
                    for featureName, count in features.items():
                        print('Feature: %s; count: %d' % (featureName, count))
            else:
                print('No significant bins.')

        if plotHashesInSubjectAndRead:
            PlotHashesInSubjectAndRead(
                query, subject, findParams=findParams,
                showInsignificant=showInsignificant, showSignificant=True,
                showBestBinOnly=showBestBinOnly, **kwargs).plotHorizontal(
                    addJitter=True)

        if showHistogram:
            plotHistogram(
                query, subject, showSignificanceCutoff=True,
                significanceFraction=findParams.significanceFraction, **kwargs)
        else:
            print('Subject %r was not matched.', subject.id)
        print()

    if showFeatures:
        # Plot landmarks and trig points horizontally.
        plotLandmarksInSequences([subject, query], **kwargs)


def scoreHeatmap(sequenceFileOrMatrix, labels, labelColors, findParams=None,
                 figureTitle=False, fileTitle=False, **kwargs):
    """
    A function to make a score heatmap.

    @param sequenceFileOrMatrix: Either a C{str} file name of a file
        containing sequences or a distance matrix as returned from
        C{light.performance.affinity}.
    @param labels: A C{list} of C{str} label names.
    @param labelColors: A C{dict} mapping each label in labels to a label
        color.
    @param findParams: A C{light.parameters.FindParameters} instance or
        C{None}, to use the default find parameters.
    @param figureTitle: If not False, a C{str} title for the figure.
    @param fileTitle: If the figure should be saved to a file, the title of the
        file where the figure is saved to.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    if isinstance(sequenceFileOrMatrix, np.ndarray):
        matrix = sequenceFileOrMatrix

    else:
        matrix = affinity.affinityMatrix(sequenceFileOrMatrix, findParams,
                                         **kwargs)

    a = plt.imshow(np.array(matrix), interpolation='nearest', cmap=plt.cm.GnBu,
                   origin='bottom')

    left, right, bottom, top = a.get_extent()

    for y, label in enumerate(labels):
        plt.text(left - 2, y, label, fontsize=15,
                 color=labelColors.get(label, 'black'))
        plt.text(left + y + 0.3, bottom - 1, label, fontsize=15,
                 rotation='vertical', color=labelColors.get(label, 'black'))
    plt.xticks([])
    plt.yticks([])
    plt.colorbar(a)
    if figureTitle:
        plt.title(figureTitle, fontsize=18)
    plt.gcf().set_size_inches(13, 10)
    if fileTitle:
        plt.savefig(fileTitle, bbox_inches='tight')
    else:
        plt.show()


def alignmentGraph(query, subject, findParams=None, createFigure=True,
                   showHistogram=True, showHorizontal=True, graphAx=None,
                   **kwargs):
    """
    Plots an alignment graph similar to the one in dark matter.

    @param query: An AARead instance of the sequence of the query.
    @param subject: An AARead instance of the sequence of the subject.
    @param findParams: An instance of C{FindParameters} or C{None}, to use the
        default find parameters.
    @param createFigure: If C{True}, create a figure and give it a title.
    @param showHistogram: If C{True}, show the histogram of the match.
    @param showHorizontal: If C{True}, show the horizontal line plot.
    @param graphAx: If not None, use this as the subplot for displaying reads.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    if createFigure:
        width = 20
        figure = plt.figure(figsize=(width, 20))

    if showHistogram:
        if showHorizontal:
            gs = gridspec.GridSpec(3, 1, height_ratios=[4, 4, 12])
            histogramAx = plt.subplot(gs[0, 0])
            horizontalAx = plt.subplot(gs[1, 0])
            graphAx = graphAx or plt.subplot(gs[2, 0])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
            histogramAx = plt.subplot(gs[0, 0])
            graphAx = graphAx or plt.subplot(gs[1, 0])
    else:
        if showHorizontal:
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
            horizontalAx = plt.subplot(gs[0, 0])
            graphAx = graphAx or plt.subplot(gs[1, 0])
        else:
            graphAx = graphAx or plt.subplot(111)

    horizontalResult = False
    findParams = findParams or FindParameters()
    if showHistogram:
        plotHistogram(query, subject, findParams=findParams,
                      readsAx=histogramAx, showMean=False, showMedian=False,
                      showStdev=False, showSignificanceCutoff=False,
                      showSignificantBins=True, **kwargs)
        histogramAx.set_title('Histogram plot (significant bins are red)',
                              fontsize=15)
    if showHorizontal:
        horizontal = PlotHashesInSubjectAndRead(
            subject, query, showSignificant=True, showInsignificant=False,
            showBestBinOnly=False, findParams=findParams, **kwargs)
        horizontalResult = horizontal.plotHorizontal(ax=horizontalAx)
        horizontalAx.set_title('Horizontal plot (top: %s, bottom: %s)' % (
                               query.id, subject.id), fontsize=15)
        horizontalAx.legend(
            handles=legendHandles(horizontalResult['namesSeen']),
            bbox_to_anchor=(0.005, 0.3, 0.5, 0.05), loc=3, ncol=2,
            borderaxespad=0.1, prop={'size': 8})

    graphAx.set_title('Bin alignments', fontsize=15)

    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    _, subjectIndex, _ = database.addSubject(subject)
    result = database.find(query, findParams, storeFullAnalysis=True)

    if subjectIndex not in result.analysis:
        graphAx.text(len(subject) / 2, 0.5, 'No match was found.',
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=15)
        maxScore = 0.0
        nSignificantBins = 0.0
        significantBins = []
    elif not result.analysis[subjectIndex]['significantBins']:
        graphAx.text(len(subject) / 2, 0.5, 'No significant bin was found.',
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=15)
        maxScore = 0.0
        nSignificantBins = 0.0
        significantBins = []
    else:
        significantBins = result.analysis[subjectIndex]['significantBins']
        maxScore = significantBins[0]['score']
        nSignificantBins = len(significantBins)

    if horizontalResult:
        cols = horizontalResult['colors']
    else:
        cols = colors.color_palette('hls', len(significantBins))

    scoreName = findParams.binScoreMethod

    totalXMin = 0
    totalXMax = len(subject)
    for i, binInfo in enumerate(significantBins):
        score = binInfo['score']
        bin_ = binInfo['bin']
        binMin = Landmark('A', len(subject), len(subject), len(subject))
        binMax = Landmark('A', 0, 0, 0)
        minFeature = None

        for match in bin_:
            lm = match['subjectLandmark']
            tp = match['subjectTrigPoint']
            # Plot the matching features.
            graphAx.plot([lm.offset, lm.offset + lm.length - 1],
                         [score, score], '-', color=COLORS[lm.symbol],
                         linewidth=3)
            graphAx.plot([tp.offset, tp.offset],
                         [score + 0.005, score - 0.005], '-',
                         color=COLORS[tp.symbol], linewidth=2)
            if lm.offset > binMax.offset:
                binMax = lm
            elif lm.offset < binMin.offset:
                binMin = lm
                minFeature = match['queryLandmark']
            if tp.offset > binMax.offset:
                binMax = tp
            elif tp.offset < binMin.offset:
                binMin = tp
                minFeature = match['queryTrigPoint']
        minBinX = binMin.offset - minFeature.offset
        maxBinX = minBinX + len(query)
        if minBinX < totalXMin:
            totalXMin = minBinX
        if maxBinX > totalXMax:
            totalXMax = maxBinX
        # plot the horizontal colored lines to indicate where the bins are
        # and the horizontal grey lines to indicate the whole query sequence.
        graphAx.plot([minBinX, maxBinX], [score, score], '-', color='grey')
        graphAx.plot([binMin.offset, binMax.offset], [score, score], '-',
                     color=cols[i], linewidth=1.5)

    # Plot vertical lines at the start and end of the subject sequence and set
    # the xlims and ylims.
    if showHorizontal:
        horizontalAx.set_xlim([totalXMin - 5, totalXMax + 5])
        horizontalAx.vlines([0.0], [0.0], [1.0], color='grey', linewidth=0.5)
        horizontalAx.vlines([len(subject)], [0.0], [1.0], color='grey',
                            linewidth=0.5)
    graphAx.vlines([0.0], [0.0], [1.0], color='grey', linewidth=0.5)
    graphAx.vlines([len(subject)], [0.0], [1.0], color='grey', linewidth=0.5)
    graphAx.set_xlim([totalXMin - 5, totalXMax + 5])
    graphAx.set_ylim([-0.01, 1.01])

    # Labels and titles
    figure.suptitle('Query: %s, Subject: %s\nmax score: %.4f, '
                    'significant bins: %d' % (
                        query.id, subject.id, maxScore, nSignificantBins),
                    fontsize=20)
    graphAx.set_ylabel(scoreName, fontsize=12)
    graphAx.set_xlabel('Sequence length (AA)', fontsize=12)
    graphAx.grid()
    figure.show()
