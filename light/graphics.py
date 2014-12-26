import matplotlib.pyplot as plt

from light.reads import ScannedRead
from light.landmarks import findLandmark
from light.trig import findTrigPoint

COLORS = {'A': 'blue',
          'B': 'cyan',
          'C': '#00BFFF',
          'S': '#819FF7',
          'P': 'red',
          'T': 'orange',
          'M': 'magenta',
          'J': 'black',
          'I': 'black',
          'O': 'black'}


def plotHistograms(lightTitlesAlignments, title):
    """
    For a lightTitlesAlignment, plot the histograms calculated by light.Result.

    @param lightTitlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param title: A C{str} sequence title that was found by light matter.
    """


def plotFeatures(read, landmarks=None, trigs=None, limitPerLandmark=None,
                 maxDistance=None, readsAx=None):
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
            maxDistance=maxDistance):
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
    readsAx.set_ylabel('Rank', fontsize=15)

    readsAx.set_xlim(0, len(read.sequence))
    readsAx.set_ylim(-0.5, count + 1)
    readsAx.grid()
