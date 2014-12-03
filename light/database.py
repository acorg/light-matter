from collections import defaultdict
import numpy as np
from scipy import stats

from light.reads import ScannedRead
from light.result import ScannedReadDatabaseResult


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

    def find(self, read):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param db: a L{light.database.ScannedRead} database in which we want to
            look up the read.
        """
        scannedRead = ScannedRead(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder(read):
                scannedRead.trigPoints.append(trigPoint)

        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            try:
                matchingKey = self.d[key]
            except KeyError:
                return
            else:
                for subjectId, subjectOffset in matchingKey.items():
                    offset = subjectOffset - landmark.offset
                    yield ScannedReadDatabaseResult(subjectId, read.id, offset)


def evaluate(found):
    """
    Evaluates whether a subject is matched significantly by a read.

    @param found: a C{list} with L{light.result.ScannedReadDatabaseResult} as
        returned by L{ScannedReadDatabase.find}.
    """
    significant = []
    for subject in found:
        for query, offsets in subject.items():
            # test significance
            # bin offsets
            hist, edges = np.histogram(offsets, bins=10)
            match = max(hist)
            t, p = stats.ttest_1samp(offsets, match)
            if p < 0.05:
                significant.append((subject, query, match))
    return significant
