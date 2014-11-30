from collections import defaultdict

from light.reads import ScannedRead


class ScannedReadDatabase(object):
    """
    Maintain a collection of reads and provide for database (index, search)
    operations on them.

    @param landmarkFinderClasses: A C{list} of landmark classes.
    @param trigPointFinderClasses: A C{list} of trig point classes.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    """
    def __init__(self, landmarkFinderClasses, trigPointFinderClasses,
                 limitPerLandmark=None, maxDistance=None):
        self.landmarkFinderClasses = landmarkFinderClasses
        self.trigPointFinderClasses = trigPointFinderClasses
        self.limitPerLandmark = limitPerLandmark
        self.maxDistance = maxDistance
        self.d = defaultdict(set)
        self.readCount = 0
        self.totalResidues = 0
        self.totalCoveredResidues = 0

        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass().find)

        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass().find)

    def addRead(self, read):
        """
        Examine a read for features and add its (landmark, trig point) pairs
        to the search dictionary.

        @param read: a C{dark.read.AARead} instance.
        """
        scannedRead = ScannedRead(read)
        self.readCount += 1
        self.totalResidues += len(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder(read):
                scannedRead.trigPoints.append(trigPoint)

        self.totalCoveredResidues += len(scannedRead.coveredIndices())

        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            self.d[key].add((read.id, landmark.offset))

    def __str__(self):
        return '%s: %d sequences, %d residues, %d hashes, %.2f%% coverage' % (
            self.__class__.__name__, self.readCount, self.totalResidues,
            len(self.d),
            float(self.totalCoveredResidues) / self.totalResidues * 100.0)
