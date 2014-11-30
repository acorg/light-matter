from collections import defaultdict

from light.reads import ScannedRead


class ScannedReadDatabase(object):
    """
    Maintain a collection of reads and provide for database (index, search)
    operations on them.

    @param landmarkFinderClasses: A C{list} of landmark classes.
    @param trigPointFinderClasses: A C{list} of trig point classes.
    """
    def __init__(self, landmarkFinderClasses, trigPointFinderClasses):
        self.landmarkFinderClasses = landmarkFinderClasses
        self.trigPointFinderClasses = trigPointFinderClasses
        self.d = defaultdict(set)

        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass().find)

        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass().find)

    def addRead(self, read):
        """
        Add (landmark, trig point) pairs to the search dictionary.

        @param read: a C{dark.read.AARead} instance.
        """
        scannedRead = ScannedRead(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder(read):
                scannedRead.trigPoints.append(trigPoint)

        for landmark, trigPoint in scannedRead.getPairs():
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            self.d[key].add((read.id, landmark.offset))
