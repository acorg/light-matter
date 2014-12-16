import sys
from os.path import basename
import matplotlib.pyplot as plt

from light.reads import ScannedRead
from light.landmarks import find as findLandmark
from light.trig import find as findTrigPoint

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


def plotFeatures(read, landmark=None, trig=None, limitPerLandmark=None,
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

    landmarks = landmark or []
    trigs = trig or []

    if len(landmarks) + len(trigs) == 0:
        print >>sys.stderr, ('You must specify either landmarks or trig '
                             'points to find.')
        sys.exit(1)

    # Make sure all landmark finders requested exist.
    landmarkFinders = []
    for landmarkFinderName in landmarks:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinders.append(landmarkFinderClass().find)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Make sure all trig point finders requested exist.
    trigFinders = []
    for trigFinderName in trigs:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinders.append(trigFinderClass().find)
        else:
            print >>sys.stderr, '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Find all landmarks and trig points on the read.
    totalCoveredResidues = 0
    scannedRead = ScannedRead(read)

    for landmarkFinder in landmarkFinders:
        for landmark in landmarkFinder(read):
            scannedRead.landmarks.append(landmark)

    for trigFinder in trigFinders:
        for trigPoint in trigFinder(read):
            scannedRead.trigPoints.append(trigPoint)

    # plot landmarks and trig point pairs.
    totalCoveredResidues += len(scannedRead.coveredIndices())
    count = 0

    for landmark, trigPoint in scannedRead.getPairs(
            limitPerLandmark=limitPerLandmark,
            maxDistance=maxDistance):
        readsAx.plot([landmark.offset, trigPoint.offset], [count, count], '-',
                     color='grey')
        landmarkColor = COLORS[landmark.hashkey()[0]]
        trigPointColor = COLORS[trigPoint.hashkey()[0]]
        readsAx.plot([landmark.offset], [count], 'o', color=landmarkColor)
        readsAx.plot([trigPoint.offset], [count], 'o', color=trigPointColor)
        count += 1

    readsAx.set_title('%s\n Length: %d, covered residues: %s' % (read.id,
                      len(read), totalCoveredResidues), fontsize=20)
    readsAx.set_ylabel('Rank', fontsize=15)

    readsAx.set_xlim(0, len(read.sequence))
    readsAx.set_ylim(-0.5, count + 1)
    readsAx.grid()
