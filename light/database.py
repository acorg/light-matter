import sys
from collections import defaultdict
from cPickle import dump, load, HIGHEST_PROTOCOL
from json import dumps
from hashlib import sha256
from operator import attrgetter


from light.reads import ScannedRead, ScannedRead as ScannedSubject
from light.result import Result


class Database(object):
    """
    Maintain a collection of sequences ("subjects") and provide for database
    (index, search) operations on them.

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
        self.subjectCount = 0
        self.totalResidues = 0
        self.totalCoveredResidues = 0
        self.subjectInfo = []

        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass())

        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass())

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
        is passed as a read instance even though in many cases it will not be
        an actual read from a sequencing run.
        """
        scannedSubject = ScannedSubject(subject)
        subjectIndex = self.subjectCount
        self.subjectCount += 1
        self.totalResidues += len(subject)
        self.subjectInfo.append((subject.id, subject.sequence))

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(subject):
                scannedSubject.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(subject):
                scannedSubject.trigPoints.append(trigPoint)

        self.totalCoveredResidues += len(scannedSubject.coveredIndices())

        for landmark, trigPoint in scannedSubject.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            self.d[key].append({"subjectIndex": subjectIndex,
                                "offset": landmark.offset,
                                })

    def __str__(self):
        return '%s: %d sequences, %d residues, %d hashes, %.2f%% coverage' % (
            self.__class__.__name__, self.subjectCount, self.totalResidues,
            len(self.d),
            float(self.totalCoveredResidues) / self.totalResidues * 100.0)

    def saveParamsAsJSON(self, fp=sys.stdout):
        """
        Save the database parameters to a file in JSON format.

        @param fp: A file pointer.
        """
        print >>fp, dumps({
            'checksum': self.checksum(),
            'landmarkFinderClasses': [
                klass.NAME for klass in self.landmarkFinderClasses],
            'trigPointFinderClasses': [
                klass.NAME for klass in self.trigPointFinderClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'subjectCount': self.subjectCount,
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

        result = Result(read, self)
        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance):
            key = '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                                landmark.offset - trigPoint.offset)
            try:
                matchingKey = self.d[key]
                landmarkLength = landmark.length
            except KeyError:
                pass
            else:
                for subjectDict in matchingKey:
                    result.addMatch(
                        {
                            'subjectOffset': subjectDict['offset'],
                            'readOffset': landmark.offset,
                        },
                        subjectDict['subjectIndex'], landmarkLength, key)
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
        @return: An instance of L{Database}.
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

    def checksum(self):
        """
        Compute a checksum for the specification and contents of this database.

        @return: A C{str} checksum.
        """
        result = sha256()

        def update(s):
            """
            Add the string representation of an object to the checksum,
            followed by a NUL.

            @param s: Anything that can be converted to a C{str} to be added to
                the checksum.
            """
            result.update(str(s) + '\0')

        # Add the (sorted) names of all landmark and trig point finder classes.
        key = attrgetter('NAME')
        for finder in sorted(self.landmarkFinders, key=key):
            update(finder.NAME)
        for finder in sorted(self.trigPointFinders, key=key):
            update(finder.NAME)

        # Add other initialization parameters.
        update(self.limitPerLandmark)
        update(self.maxDistance)

        # Add all subject info.
        for subjectId, subjectSequence in self.subjectInfo:
            update(subjectId)
            update(subjectSequence)

        # Add all hash keys and their content. We do this in sorted order
        # for the keys of self.d and then within each dict in the self.d
        # values we make sure we have the sorted order of those keys
        # too. This makes sure the checksum is invariant with respect to
        # the order in which Python enumerates its dictionary keys.
        for key in sorted(self.d.keys()):
            update(key)
            for details in self.d[key]:
                for detailKey in sorted(details.keys()):
                    update(detailKey)
                    update(details[detailKey])

        return result.hexdigest()
