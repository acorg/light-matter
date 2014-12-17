import sys
from collections import defaultdict
from cPickle import dump, load, HIGHEST_PROTOCOL
from json import dumps


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
        self.d = defaultdict(list)
        self.readCount = 0
        self.totalResidues = 0
        self.totalCoveredResidues = 0
        self.readInfo = []

        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass())

        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass())

    def addRead(self, read):
        """
        Examine a read for features and add its (landmark, trig point) pairs
        to the search dictionary.

        @param read: a C{dark.read.AARead} instance.
        """
        scannedRead = ScannedRead(read)
        readIndex = self.readCount
        self.readCount += 1
        self.totalResidues += len(read)
        self.readInfo.append((read.id, len(read)))

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(read):
                scannedRead.trigPoints.append(trigPoint)

        self.totalCoveredResidues += len(scannedRead.coveredIndices())

        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            self.d[key].append({"subjectIndex": readIndex,
                                "offset": landmark.offset,
                                "length": len(read),
                                })

    def __str__(self):
        return '%s: %d sequences, %d residues, %d hashes, %.2f%% coverage' % (
            self.__class__.__name__, self.readCount, self.totalResidues,
            len(self.d),
            float(self.totalCoveredResidues) / self.totalResidues * 100.0)

    def saveParamsAsJSON(self, fp=sys.stdout):
        """
        Save the database parameters to a file in JSON format.

        @param fp: A file pointer.
        """
        print >>fp, dumps({
            'landmarkFinderClasses': [
                klass.NAME for klass in self.landmarkFinderClasses],
            'trigPointFinderClasses': [
                klass.NAME for klass in self.trigPointFinderClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'readCount': self.readCount,
            'totalResidues': self.totalResidues,
            'totalCoveredResidues': self.totalCoveredResidues,
        })

    def find(self, read):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: a C{dark.read.AARead} instance.
        """
        scannedRead = ScannedRead(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(read):
                scannedRead.trigPoints.append(trigPoint)

        result = ScannedReadDatabaseResult(read, self)
        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            try:
                matchingKey = self.d[key]
            except KeyError:
                pass
            else:
                for subjectDict in matchingKey:
                    result.addMatch(
                        {
                            'subjectOffset': subjectDict['offset'],
                            'readOffset': landmark.offset,
                        },
                        subjectDict['subjectIndex'])
        result.finalize()
        return result

    def save(self, fp=sys.stdout):
        """
        Save the database to a file.

        @param fp: A file pointer.
        """
        dump(self, fp, protocol=HIGHEST_PROTOCOL)

    @staticmethod
    def load(fp=sys.stdin):
        """
        Load a database from a file.

        @param fp: A file pointer.
        @return: An instance of L{ScannedReadDatabase}.
        """
        # NOTE: We're using pickle, which isn't considered secure. But running
        # other people's Python code in general shouldn't be considered
        # secure, either. Make sure you don't load saved databases from
        # untrusted sources. We could write this using JSON, but the set
        # objects in self.d are not serializable and I don't want to convert
        # each to a tuple due to concerns about memory usage when databases
        # grow large. Anyway, for now let's proceed with pickle and caution.
        # See google for more on Python and pickle security.
        return load(fp)
