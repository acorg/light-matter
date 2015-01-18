from random import uniform
import matplotlib.pyplot as plt
from matplotlib import patches
from operator import attrgetter

from light.reads import ScannedRead
from light.database import Database
from light.trig import findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.landmarks import findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.features import Landmark

from dark.dimension import dimensionalIterator

ALL_FEATURES = [
    (feature.SYMBOL, feature.NAME) for feature in
    sorted(ALL_LANDMARK_FINDER_CLASSES | ALL_TRIG_FINDER_CLASSES,
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


def plotHistogram(query, subject, landmarks=None, trigPoints=None,
                  limitPerLandmark=None, maxDistance=None, minDistance=None,
                  aboveMeanTreshold=None, bucketFactor=1, readsAx=None):
    """
    A function which plots a histogram of matching hash distances.

    @param query: an AARead instance of the sequence of the query.
    @param subject: an AARead instance of the sequence of the subject.
    @param landmarks: a C{list} of C{str} of landmark finder names.
    @param trigPoints: a C{list} of C{str} of trig finder names.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param aboveMeanThreshold: A numeric amount by which the maximum delta
        count in a bucket must exceed the mean bucket count for that
        maximum bucket count to be considered significant.
    @para bucketFactor: A C{int} factor by which the distance between landmark
        and trig point is divided.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    """
    fig = plt.figure()
    readsAx = readsAx or fig.add_subplot(111)

    landmarks = landmarks or []
    trigs = trigPoints or []

    if len(landmarks) + len(trigs) == 0:
        raise ValueError('You must specify either landmarks or trig points to '
                         'find.')

    landmarkFinderClasses = []
    for landmarkFinderName in landmarks:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinderClasses.append(landmarkFinderClass)
        else:
            print 'Could not find landmark finder %r.' % (
                landmarkFinderName)

    # Make sure all trig point finders requested exist.
    trigFinderClasses = []
    for trigFinderName in trigs:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinderClasses.append(trigFinderClass)
        else:
            print 'Could not find trig point finder %r.' % (
                trigFinderName)

    database = Database(landmarkFinderClasses, trigFinderClasses,
                        limitPerLandmark, maxDistance, minDistance,
                        bucketFactor)
    database.addSubject(subject)

    result = database.find(query, aboveMeanTreshold, storeAnalysis=True)
    hist = result.analysis[0]['histogram']
    bins = result.analysis[0]['histogramBuckets']

    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    readsAx.bar(center, hist, align='center', width=width)

    readsAx.set_title('%s against %s' % (query.id, subject.id))

    readsAx.set_xlabel("Offsets (database-read)")

    readsAx.xaxis.tick_bottom()


def plotFeatures(read, landmarks=None, trigs=None, limitPerLandmark=None,
                 maxDistance=None, minDistance=None, readsAx=None):
    """
    A function which plots the positions of landmark and trigpoint pairs on a
    sequence.

    @param read: A C{dark.reads.Read} instance.
    @param landmark: a C{list} of C{str} of landmark finder names.
    @param trig: a C{list} of C{str} of trig finder names.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    """
    width = 20
    fig = plt.figure(figsize=(width, 20))
    readsAx = readsAx or fig.add_subplot(111)

    landmarks = landmarks or []
    trigs = trigs or []

    if len(landmarks) + len(trigs) == 0:
        raise ValueError('You must specify either landmarks or trig points to '
                         'find.')

    # Make sure all landmark finders requested exist.
    landmarkFinders = []
    for landmarkFinderName in landmarks:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinders.append(landmarkFinderClass().find)
        else:
            raise ValueError('Could not find landmark finder %r.' % (
                             landmarkFinderName))

    # Make sure all trig point finders requested exist.
    trigFinders = []
    for trigFinderName in trigs:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinders.append(trigFinderClass().find)
        else:
            raise ValueError('Could not find trig point finder %r.' % (
                             trigFinderName))

    # Find all landmarks and trig points on the read.
    scannedRead = ScannedRead(read)

    for landmarkFinder in landmarkFinders:
        for landmark in landmarkFinder(read):
            scannedRead.landmarks.append(landmark)

    for trigFinder in trigFinders:
        for trigPoint in trigFinder(read):
            scannedRead.trigPoints.append(trigPoint)

    # plot landmarks and trig point pairs.
    totalCoveredResidues = len(scannedRead.coveredIndices())
    count = 0

    for landmark, trigPoint in scannedRead.getPairs(
            limitPerLandmark=limitPerLandmark,
            maxDistance=maxDistance, minDistance=minDistance):
        readsAx.plot([landmark.offset, trigPoint.offset], [count, count], '-',
                     color='grey')
        landmarkColor = COLORS[landmark.symbol]
        trigPointColor = COLORS[trigPoint.symbol]
        readsAx.plot([landmark.offset, landmark.offset + landmark.length],
                     [count, count], '-', color=landmarkColor, linewidth=4)
        readsAx.plot([trigPoint.offset, trigPoint.offset + trigPoint.length],
                     [count, count], '-', color=trigPointColor, linewidth=4)
        count += 1

    readsAx.set_title('%s\n Length: %d, covered residues: %s' % (read.id,
                      len(read), totalCoveredResidues), fontsize=20)
    readsAx.set_ylabel('Rank', fontsize=10)

    readsAx.set_xlim(0, len(read.sequence))
    readsAx.set_ylim(-0.5, count + 1)
    readsAx.grid()

    return scannedRead, count


def plotFeaturePanel(reads, landmarks=None, trigs=None, limitPerLandmark=None,
                     maxDistance=None, minDistance=None):
    """
    Plot a panel of feature plots from plotFeatures.

    @param reads: A C{dark.fasta.FastaReads} instance.
    @param landmark: a C{list} of C{str} of landmark finder names.
    @param trig: a C{list} of C{str} of trig finder names.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    """
    cols = 5
    rows = 40
    # rows = int(len(reads) / cols) + (0 if len(reads) % cols == 0 else 1)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    maxY = 0
    maxX = 0

    coords = dimensionalIterator((rows, cols))

    for i, read in enumerate(reads):
        row, col = coords.next()
        print '%d: %s' % (i, read.id)

        graphInfo, y = plotFeatures(read, landmarks, trigs,
                                    limitPerLandmark, maxDistance,
                                    minDistance, readsAx=ax[row][col])

        totalCoveredResidues = len(graphInfo.coveredIndices())
        plotTitle = ('%s\n Length: %d, covered residues: %s' % (read.id,
                     len(read), totalCoveredResidues))

        ax[row][col].set_title(plotTitle, fontsize=10)
        if y > maxY:
            maxY = y
        if len(read) > maxX:
            maxX = len(read)

    # Post-process graphs to adjust axes, etc.
    coords = dimensionalIterator((rows, cols))
    for read in reads:
        row, col = coords.next()
        a = ax[row][col]
        # a.set_ylim([-0.5, maxY + 0.5])
        a.set_yticks([])
        a.set_xticks([])
        a.set_ylabel('')

    # Hide the final panel graphs (if any) that have no content. We do this
    # because the panel is a rectangular grid and some of the plots at the
    # end of the last row may be unused.
    for row, col in coords:
        ax[row][col].axis('off')

    plt.subplots_adjust(hspace=0.4)
    figure.suptitle('X: 0 to %d, Y: 0 to %d (rank)' % (maxX, maxY),
                    fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)
    figure.show()


def plotFeatureSquare(read, landmarks=None, trigs=None, limitPerLandmark=None,
                      maxDistance=None, minDistance=None, readsAx=None):
    """
    Plot the positions of landmark and trigpoint pairs on a sequence in a
    square.

    @param read: A C{dark.reads.Read} instance.
    @param landmark: a C{list} of C{str} of landmark finder names.
    @param trig: a C{list} of C{str} of trig finder names.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    """
    fig = plt.figure(figsize=(15, 15))
    readsAx = readsAx or fig.add_subplot(111)
    landmarks = landmarks or []
    trigs = trigs or []

    if len(landmarks) + len(trigs) == 0:
        raise ValueError('You must specify either landmarks or trig points to '
                         'find.')

    # Make sure all landmark finders requested exist.
    landmarkFinders = []
    for landmarkFinderName in landmarks:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinders.append(landmarkFinderClass().find)
        else:
            raise ValueError('Could not find landmark finder %r.' % (
                             landmarkFinderName))

    # Make sure all trig point finders requested exist.
    trigFinders = []
    for trigFinderName in trigs:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinders.append(trigFinderClass().find)
        else:
            raise ValueError('Could not find trig point finder %r.' % (
                             trigFinderName))

    # Find all landmarks and trig points on the read.
    scannedRead = ScannedRead(read)

    for trigFinder in trigFinders:
        for trigPoint in trigFinder(read):
            scannedRead.trigPoints.append(trigPoint)
    for landmarkFinder in landmarkFinders:
        for landmark in landmarkFinder(read):
            scannedRead.landmarks.append(landmark)

    # Plot a light grey diagonal line, bottom left to top right.
    readsAx.plot([0, len(read.sequence)], [0, len(read.sequence)], '-',
                 color='#f0f0f0', linewidth=1)

    scatterX = []
    scatterY = []
    scatterColors = []
    namesSeen = set()
    landmarks = set()

    for landmark, trigPoint in scannedRead.getPairs(
            limitPerLandmark=limitPerLandmark,
            maxDistance=maxDistance, minDistance=minDistance):
        # Add jitter to the Y offset so we can see more trig points that
        # occur at the same offset.
        scatterX.append(landmark.offset)
        scatterY.append(trigPoint.offset + uniform(-0.4, 0.4))
        if isinstance(trigPoint, Landmark):
            # Color the trig points that are actually landmarks a light grey.
            scatterColors.append((0.9, 0.9, 0.9))
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
    totalCoveredResidues = len(scannedRead.coveredIndices())
    readsAx.set_title('%s\n Length: %d, covered residues: %s' % (read.id,
                      len(read), totalCoveredResidues), fontsize=20)
    readsAx.set_xlabel('Offset (landmarks)', fontsize=16)
    readsAx.set_ylabel('Offset (trig points)', fontsize=16)
    readsAx.set_xlim(-0.2, len(read.sequence) + 0.2)
    readsAx.set_ylim(0, len(read.sequence))
    readsAx.legend(handles=legendHandles(namesSeen), loc=2)
    readsAx.grid()

    return scannedRead
