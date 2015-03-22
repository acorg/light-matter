from random import uniform
import matplotlib.pyplot as plt
from matplotlib import patches
from operator import attrgetter
from textwrap import fill
from collections import defaultdict

from dark.dimension import dimensionalIterator

from light.database import DatabaseSpecifier
from light.trig import ALL_TRIG_CLASSES
from light.landmarks import ALL_LANDMARK_CLASSES
from light.features import Landmark
from light.colors import colors

from dark.fasta import FastaReads
from dark.reads import AARead

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
                      equalizeYAxes=False):
    """
    Create a rectangular panel of plots.

    @param rows: The C{int} number of rows that will appear in the panel.
    @param cols: The C{int} number of columns that will appear in the panel.
    @param title: A C{str} title for the figure.
    @param makeSubPlot: A function that accepts two arguments (an C{int} row
        and column) and which creates a sub-plot and returns a C{dict}
        containing information about the created plot or C{None} if no plot
        should be shown in the given (row, col) location.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param equalizeYAxes: if C{True}, adjust the Y axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    """
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    subplots = {}

    for row, col in dimensionalIterator((rows, cols)):
        subplots[(row, col)] = makeSubPlot(row, col, ax[row][col])

    if equalizeXAxes or equalizeYAxes:
        nonEmpty = filter(lambda x: x, subplots.itervalues())
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
    for (row, col), subplot in subplots.iteritems():
        a = ax[row][col]
        if subplot:
            try:
                subTitle = subplots[(row, col)]['title']
            except KeyError:
                # No title, no problem.
                pass
            else:
                a.set_title(subTitle, fontsize=10)
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
    figure.set_size_inches(15, 3 * rows, forward=True)
    figure.show()


def plotHistogramPanel(sequences, equalizeXAxes=True, equalizeYAxes=False,
                       significanceFraction=None, **kwargs):
    """
    Plot a panel of histograms of matching hash offset deltas between
    all pairs of passed sequences.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @param equalizeXAxes: if C{True}, adjust the X axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @param equalizeYAxes: if C{True}, adjust the Y axis on each sub-plot
        to cover the same range (the maximum range of all sub-plots).
    @return: The C{light.result.Result} from running the database find.
    """
    if isinstance(sequences, basestring):
        reads = list(FastaReads(sequences, readClass=AARead))
    else:
        reads = list(sequences)
    nReads = len(reads)
    database = DatabaseSpecifier().getDatabaseFromKeywords(subjects=reads,
                                                           **kwargs)

    def makeSubPlot(row, col, ax):
        """
        @param row: The C{int} panel row.
        @param col: The C{int} panel column.
        @param ax: The matplotlib axis for the sub-plot.
        """
        if row <= col:
            return plotHistogram(sequences[row], sequences[col],
                                 significanceFraction=significanceFraction,
                                 readsAx=ax, database=database)

    return _rectangularPanel(
        nReads, nReads, 'Histogram panel', makeSubPlot,
        equalizeXAxes=equalizeXAxes, equalizeYAxes=equalizeYAxes)


def plotHistogram(query, subject, significanceFraction=None, readsAx=None,
                  **kwargs):
    """
    Plot a histogram of matching hash offset deltas between a query and a
    subject.

    @param query: an AARead instance of the sequence of the query.
    @param subject: an AARead instance of the sequence of the subject.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    @return: The C{light.result.Result} from running the database find.
    """
    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    subjectIndex = database.addSubject(subject)
    result = database.find(query, significanceFraction, storeFullAnalysis=True)

    try:
        histogram = result.analysis[subjectIndex]['histogram']
    except KeyError:
        print 'Query %r and subject %r had no hashes in common.' % (
            query.id, subject.id)
    else:
        counts = [len(bin) for bin in histogram.bins]
        nBins = len(histogram.bins)
        width = (histogram.max - histogram.min) / float(nBins)
        center = [histogram.min + (i + 0.5) * width for i in xrange(nBins)]
        title = fill('%s vs %s' % (query.id, subject.id), FILL_WIDTH)

        if readsAx is None:
            fig = plt.figure()
            readsAx = fig.add_subplot(111)
            readsAx.set_title(title)
            readsAx.set_xlabel('Offsets (database-read)')
            readsAx.xaxis.tick_bottom()

        readsAx.bar(center, counts, align='center', width=width)
        return {
            'minX': histogram.min,
            'maxX': histogram.max,
            'minY': 0,
            'maxY': max(counts),
            'result': result,
            'title': title,
        }


def plotFeatureSquare(read, significanceFraction=None, readsAx=None, **kwargs):
    """
    Plot the positions of landmark and trigpoint pairs on a sequence in a
    square.

    @param read: A C{dark.reads.Read} instance.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    fig = plt.figure(figsize=(15, 15))
    readsAx = readsAx or fig.add_subplot(111)

    database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    result = database.find(read, significanceFraction, storeFullAnalysis=True)
    scannedQuery = result.scannedQuery

    # Plot a light grey diagonal line, bottom left to top right.
    readsAx.plot([0, len(read.sequence)], [0, len(read.sequence)], '-',
                 color='#f0f0f0', linewidth=1)

    scatterX = []
    scatterY = []
    scatterColors = []
    namesSeen = set()
    landmarks = set()

    for landmark, trigPoint in database.getScannedPairs(scannedQuery):
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
        readsAx.plot([offset, offset + length], [offset, offset], '-',
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
    readsAx.legend(handles=legendHandles(namesSeen), loc=2)
    readsAx.grid()

    return scannedQuery


class PlotHashesInSubjectAndRead(object):
    """
    A class which plots visualisations of the hashes in subject and query.
    It collects three types of hashes: 1) Hashes in the query that don't
    match in the subject. 2) Hashes in the subject that don't match in the
    query. 3) Hashes that match in subject and query. These can subsequently
    be plotted.

    @param query: an AARead instance of the sequence of the query.
    @param subject: an AARead instance of the sequence of the subject.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param kwargs: See C{database.DatabaseSpecifier.getDatabaseFromKeywords}
        for additional keywords, all of which are optional.
    """
    def __init__(self, query, subject, significanceFraction=None, **kwargs):
        self.query = query
        self.subject = subject

        database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        subjectIndex = database.addSubject(subject)
        self.result = database.find(self.query, significanceFraction,
                                    storeFullAnalysis=True)
        if subjectIndex in self.result.analysis:
            analysis = self.result.analysis[subjectIndex]
            self.score = analysis['bestScore']
            self.significantBinCount = len(analysis['significantBins'])
            self.bins = analysis['histogram'].bins
        else:
            self.score = 0.0
            self.significantBinCount = 0
            self.bins = {}

        self.queryHashes = self.result.nonMatchingHashes
        self.subjectHashes = database.getHashes(database.scan(subject))
        self.matchingHashes = defaultdict(list)

        for bin_ in self.bins:
            for match in bin_:
                hash_ = database.hash(match['landmark'], match['trigPoint'])
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

        for hashInfo in self.queryHashes.itervalues():
            for offset in hashInfo['landmarkOffsets']:
                readsAx.plot(offset + uniform(-0.4, 0.4), 0, 'o',
                             markerfacecolor='black', markeredgecolor='white')

        for hashInfo in self.subjectHashes.itervalues():
            for offset in hashInfo['landmarkOffsets']:
                readsAx.plot(0, offset + uniform(-0.4, 0.4), 'o',
                             markerfacecolor='black', markeredgecolor='white')

        nonEmptyBins = filter(lambda b: len(b), self.bins)
        binColors = colors.color_palette('hls', len(nonEmptyBins))
        for bin_, binColor in zip(nonEmptyBins, binColors):
            for match in bin_:
                readsAx.plot(
                    match['subjectLandmarkOffset'] + uniform(-0.4, 0.4),
                    match['queryLandmarkOffset'] + uniform(-0.4, 0.4),
                    'o', markerfacecolor=binColor, markeredgecolor='white')

        firstTitleLine = fill('Hashes from matching %s against %s' % (
            self.query.id, self.subject.id), FILL_WIDTH)
        readsAx.set_title(
            '%s\nScore: %.4f, HSPs: %d\n'
            'Hashes: matching=%d, subject-only=%d, query-only=%d' % (
                firstTitleLine, self.score or 0.0, self.significantBinCount,
                len(self.matchingHashes), len(self.subjectHashes),
                len(self.queryHashes)))
        readsAx.set_ylabel('Query: %s' % self.query.id)
        readsAx.set_xlabel('Subject: %s' % self.subject.id)
        readsAx.set_xlim(-0.5, len(self.subject))
        readsAx.set_ylim(-0.5, len(self.query))
        readsAx.grid()

    def plotHorizontal(self, ax=None):
        """
        Plot a graph with the subject and query drawn as two horizontal lines.

        @param ax: If not C{None}, use this as the subplot.
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
        # handling, and landmarks are handled befor trig points.

        fig = plt.figure(figsize=(15, 10))
        ax = ax or fig.add_subplot(111)
        namesSeen = set()
        qyY = 0.0
        sjY = 1.0
        tpY = 0.015
        missY = 0.1
        maxLen = max([len(self.query), len(self.subject)])
        verticalPad = 0.1
        horizontalPad = 0.02 * maxLen

        # A line for query hashes that match the subject.
        ax.plot([0.0, len(self.query)], [qyY, qyY],
                '-', color='#bbbbbb')

        # A lower line for query hashes that don't match the subject.
        ax.plot([0.0, len(self.query)], [qyY - missY, qyY - missY],
                '-', color='#cccccc')

        # A line for subject hashes that match the query.
        ax.plot([0.0, len(self.subject)], [sjY, sjY],
                '-', color='#bbbbbb')

        # A higher line for subject hashes that don't match the query.
        ax.plot([0.0, len(self.subject)], [sjY + missY, sjY + missY],
                '-', color='#cccccc')

        # Keep track of features plotted in the query and subject, so we
        # don't plot them in the missed feature areas. Similar for the
        # diagonal lines showing matches.
        qyPlotted = set()
        sjPlotted = set()
        diagPlotted = set()

        # Plot matching hashes on the query and subject horizontal lines,
        # as well as diagonal lines between the two (in colors that
        # correspond to histogram bins).

        nonEmptyBins = filter(lambda b: len(b), self.bins)
        binColors = colors.color_palette('hls', len(nonEmptyBins))

        for bin_, binColor in zip(nonEmptyBins, binColors):
            for match in bin_:
                lm = match['landmark']
                tp = match['trigPoint']
                lmName = lm.name
                tpName = tp.name
                lmColor = COLORS[lm.symbol]
                tpColor = COLORS[tp.symbol]
                qyLmOffset = match['queryLandmarkOffset']
                qyTpOffset = match['queryTrigPointOffset']
                sjLmOffset = match['subjectLandmarkOffset']
                namesSeen.update([lmName, tpName])

                # Landmark in the query.
                key = (lmName, qyLmOffset, lm.length)
                if key not in qyPlotted:
                    qyPlotted.add(key)
                    plt.plot([qyLmOffset, qyLmOffset + lm.length], [qyY, qyY],
                             '-', color=lmColor, linewidth=3)

                # Trig point in the query.
                key = (tpName, qyTpOffset)
                if key not in qyPlotted:
                    qyPlotted.add(key)
                    plt.plot([qyTpOffset, qyTpOffset], [qyY + tpY, qyY - tpY],
                             '-', color=tpColor, linewidth=3)

                # Landmark in the subject.
                key = (lmName, sjLmOffset, lm.length)
                if key not in sjPlotted:
                    sjPlotted.add(key)
                    plt.plot([sjLmOffset, sjLmOffset + lm.length], [sjY, sjY],
                             '-', color=lmColor, linewidth=3)

                # Trig point in the subject.
                #
                # Note that we don't know the subject trig point offset for
                # sure. We assume here that the distance between landmark
                # and trig point in the subject is the same as it is in the
                # query. That could be slightly wrong, depending on the
                # landmark to trig point distance scaling in database.hash().
                sjTpOffset = sjLmOffset + (qyTpOffset - qyLmOffset)
                key = (tpName, sjTpOffset)
                if key not in sjPlotted:
                    sjPlotted.add(key)
                    plt.plot([sjTpOffset, sjTpOffset], [sjY + tpY, sjY - tpY],
                             '-', color=lmColor, linewidth=3)

                # Add x-axis offset jitter to diagonal line plotting, so we
                # can see more of them in case of overlap.
                xJitter = uniform(-0.2, 0.2)

                # Dashed diagonal line connecting trig point in query and
                # subject.  Plot this before the diagonal landmark lines as
                # it's arguably less important (there's a chance we will
                # paint over it with the landmark lines).
                key = (tpName, qyTpOffset, sjTpOffset)
                if key not in diagPlotted:
                    diagPlotted.add(key)
                    plt.plot([qyTpOffset + xJitter, sjTpOffset + xJitter],
                             [qyY, sjY], '--', color=binColor)

                # Solid diagonal line connecting start of landmark in query
                # and subject.
                key = (lmName, qyLmOffset, sjLmOffset)
                if key not in diagPlotted:
                    diagPlotted.add(key)
                    plt.plot([qyLmOffset + xJitter, sjLmOffset + xJitter],
                             [qyY, sjY], '-', color=binColor)

                # Don't plot the end of the landmark, to reduce clutter.

        # Query-only hashes, plotted just below (-missY) the query line.
        for hashInfo in self.queryHashes.itervalues():
            lm = hashInfo['landmark']
            tp = hashInfo['trigPoint']
            lmName = lm.name
            tpName = tp.name
            lmColor = COLORS[lm.symbol]
            tpColor = COLORS[tp.symbol]
            namesSeen.update([lmName, tpName])
            for lmOffset, tpOffset in zip(
                    hashInfo['landmarkOffsets'], hashInfo['trigPointOffsets']):
                # Landmark.
                key = (lmName, lmOffset, lm.length)
                if key not in qyPlotted:
                    qyPlotted.add(key)
                    plt.plot([lmOffset, lmOffset + lm.length],
                             [qyY - missY, qyY - missY],
                             '-', color=lmColor, linewidth=3)
                # Trig point.
                key = (tpName, tpOffset)
                if key not in qyPlotted:
                    qyPlotted.add(key)
                    plt.plot([tpOffset, tpOffset],
                             [qyY - missY + tpY, qyY - missY - tpY],
                             '-', color=tpColor, linewidth=3)

        # Subject-only hashes, plotted just above (+missY) the subject line.
        for hashInfo in self.subjectHashes.itervalues():
            lm = hashInfo['landmark']
            tp = hashInfo['trigPoint']
            lmName = lm.name
            tpName = tp.name
            lmColor = COLORS[lm.symbol]
            tpColor = COLORS[tp.symbol]
            namesSeen.update([lmName, tpName])
            for lmOffset, tpOffset in zip(
                    hashInfo['landmarkOffsets'], hashInfo['trigPointOffsets']):
                # Landmark.
                key = (lmName, lmOffset, lm.length)
                if key not in sjPlotted:
                    sjPlotted.add(key)
                    plt.plot([lmOffset, lmOffset + lm.length],
                             [sjY + missY, sjY + missY],
                             '-', color=lmColor, linewidth=3)
                # Trig point.
                key = (tpName, tpOffset)
                if key not in sjPlotted:
                    sjPlotted.add(key)
                    plt.plot([tpOffset, tpOffset],
                             [sjY + missY + tpY, sjY + missY - tpY],
                             '-', color=tpColor, linewidth=3)

        ax.legend(handles=legendHandles(namesSeen),
                  bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3, ncol=2,
                  borderaxespad=0.5)
        ax.set_xlabel(fill('Subject (top): %s\n\n\nQuery (bottom): %s' %
                           (self.subject.id, self.query.id)), fontsize=17)
        ax.set_ylim(qyY - missY - verticalPad, sjY + missY + verticalPad)
        ax.set_yticks([])
        ax.set_xlim(-horizontalPad, maxLen + horizontalPad)


def plotLandmarksInSequences(sequences, **kwargs):
    """
    Plot the positions of landmarks and trig points on many sequences, with
    sequences stacked above each other.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    if isinstance(sequences, basestring):
        reads = list(FastaReads(sequences, readClass=AARead))
    else:
        reads = list(sequences)

    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    fig = plt.figure(figsize=(15, len(reads) / 3))
    ax = fig.add_subplot(111)
    namesSeen = set()
    maxLen = 0

    for i, read in enumerate(reads):
        readLen = len(read)
        if readLen > maxLen:
            maxLen = readLen
        plt.plot([0, len(read.sequence)], [i, i], '-', linewidth=0.5,
                 color='grey')
        scannedRead = db.scan(read)
        # Landmarks are drawn as colored horizontal lines.
        for landmark in scannedRead.landmarks:
            namesSeen.add(landmark.name)
            plt.plot([landmark.offset, landmark.offset + landmark.length],
                     [i, i], '-', color=COLORS[landmark.symbol], linewidth=2)
        # Trig points are drawn as small colored vertical lines.
        for trigPoint in scannedRead.trigPoints:
            namesSeen.add(trigPoint.name)
            plt.plot([trigPoint.offset, trigPoint.offset],
                     [i - 0.125, i + 0.125], '-',
                     color=COLORS[trigPoint.symbol], linewidth=2)

    ax.set_title('Landmarks and features')
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.xaxis.grid()
    ax.set_ylim(-0.1, i + 0.1)
    ax.set_xlim(0, maxLen)
    # Add a legend above left on the plot.
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
    colorbar = plt.colorbar(ticks=[range(0,
                            max(confusionMatrix.max(axis=1)) + 1)])
    colorbar.outline.remove()
    plt.tick_params(axis='x', which='both', bottom='off', top='off',
                    labelbottom='off', labeltop='on')
    plt.tick_params(axis='y', which='both', left='off', right='off',
                    labelbottom='on')
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
