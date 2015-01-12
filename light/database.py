import sys
from collections import defaultdict
from ujson import dump, dumps, load
from hashlib import sha256
from operator import attrgetter


from light.reads import ScannedRead, ScannedRead as ScannedSubject
from light.result import Result
from light.landmarks import findLandmark
from light.trig import findTrigPoint


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
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param bucketFactor: A C{int} factor by which the distance between
        landmark and trig point is divided, to influence sensitivity.
    """

    # The default amount by which the maximum delta count in a bucket must
    # exceed the mean bucket count for that maximum bucket count to be
    # considered significant.
    ABOVE_MEAN_THRESHOLD_DEFAULT = 15

    def __init__(self, landmarkFinderClasses, trigPointFinderClasses,
                 limitPerLandmark=None, maxDistance=None, minDistance=None,
                 bucketFactor=1):
        self.landmarkFinderClasses = landmarkFinderClasses
        self.trigPointFinderClasses = trigPointFinderClasses
        self.limitPerLandmark = limitPerLandmark
        self.maxDistance = maxDistance
        self.minDistance = minDistance
        # The factor by which the distance of landmark and trigpoint pairs is
        # divided, to influence sensitivity.
        self.bucketFactor = bucketFactor
        # It may look like self.d should be a defaultdict(list). But that
        # will not work because a database JSON save followed by a load
        # will restore the defaultdict as a vanilla dict.
        self.d = {}
        self.subjectCount = 0
        self.totalResidues = 0
        self.totalCoveredResidues = 0
        self.subjectInfo = []
        self._checksum = None
        # Create instances of the landmark and trig point finder classes.
        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass())
        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass())

    def key(self, landmark, trigPoint):
        """
        Compute a key to store information about a landmark / trig point
        association for a read.

        @param landmark: A C{light.features.Landmark} instance.
        @param trigPoint: A C{light.features.TrigPoint} instance.
        @return: A C{str} key based on the landmark, the trig point,
            and the distance between them.
        """
        distance = ((landmark.offset - trigPoint.offset)
                    // self.bucketFactor)
        return '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                             distance)

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        """
        # Invalidate the stored checksum (if any).
        self._checksum = None
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
                maxDistance=self.maxDistance, minDistance=self.minDistance):
            key = self.key(landmark, trigPoint)
            value = {
                'subjectIndex': subjectIndex,
                'offset': landmark.offset,
            }
            try:
                self.d[key].append(value)
            except KeyError:
                self.d[key] = [value]

    def __str__(self):
        return '%s: %d sequences, %d residues, %d hashes, %.2f%% coverage' % (
            self.__class__.__name__, self.subjectCount, self.totalResidues,
            len(self.d),
            float(self.totalCoveredResidues) / self.totalResidues * 100.0)

    def saveParamsAsJSON(self, fp=sys.stdout):
        """
        Save the database parameters to a file in JSON format.

        @param fp: A file pointer.
        @return: The C{fp} we were passed (this is useful in testing).
        """
        print >>fp, dumps({
            'checksum': self.checksum(),
            'landmarkFinderClasses': [
                klass.NAME for klass in self.landmarkFinderClasses],
            'trigPointFinderClasses': [
                klass.NAME for klass in self.trigPointFinderClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'minDistance': self.minDistance,
            'subjectCount': self.subjectCount,
            'totalResidues': self.totalResidues,
            'totalCoveredResidues': self.totalCoveredResidues,
            'bucketFactor': self.bucketFactor,
        })

        return fp

    def find(self, read, aboveMeanThreshold=None, storeAnalysis=False):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: a C{dark.read.AARead} instance.
        @param aboveMeanThreshold: A numeric amount by which the maximum delta
            count in a bucket must exceed the mean bucket count for that
            maximum bucket count to be considered significant.
        @param storeAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if aboveMeanThreshold is None:
            aboveMeanThreshold = self.ABOVE_MEAN_THRESHOLD_DEFAULT

        scannedRead = ScannedRead(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(read):
                scannedRead.trigPoints.append(trigPoint)

        matches = defaultdict(list)

        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance, minDistance=self.minDistance):
            key = self.key(landmark, trigPoint)
            try:
                subjects = self.d[key]
            except KeyError:
                pass
            else:
                for subject in subjects:
                    matches[subject['subjectIndex']].append({
                        'distance': trigPoint.offset - landmark.offset,
                        'landmarkLength': landmark.length,
                        'landmarkName': landmark.name,
                        'readOffset': landmark.offset,
                        'subjectOffset': subject['offset'],
                        'trigPointName': trigPoint.name,
                    })

        return Result(read, matches, aboveMeanThreshold,
                      storeAnalysis=storeAnalysis)

    def save(self, fp=sys.stdout):
        """
        Save the database to a file.

        @param fp: A file pointer.
        """
        state = {
            '_checksum': self.checksum(),
            'landmarkFinderClassNames': [klass.NAME for klass in
                                         self.landmarkFinderClasses],
            'trigPointFinderClassNames': [klass.NAME for klass in
                                          self.trigPointFinderClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'minDistance': self.minDistance,
            'd': self.d,
            'subjectCount': self.subjectCount,
            'totalResidues': self.totalResidues,
            'totalCoveredResidues': self.totalCoveredResidues,
            'subjectInfo': self.subjectInfo,
            'bucketFactor': self.bucketFactor
        }
        dump(state, fp)

    @staticmethod
    def load(fp=sys.stdin):
        """
        Load a database from a file.

        @param fp: A file pointer.
        @return: An instance of L{Database}.
        """
        state = load(fp)

        landmarkFinderClasses = []
        for landmarkClassName in state['landmarkFinderClassNames']:
            klass = findLandmark(landmarkClassName)
            if klass:
                landmarkFinderClasses.append(klass)
            else:
                print >>sys.stderr, (
                    'Could not find landscape finder class %r! Has that '
                    'class been renamed or removed?' % landmarkClassName)
                sys.exit(1)

        trigPointFinderClasses = []
        for trigPointClassName in state['trigPointFinderClassNames']:
            klass = findTrigPoint(trigPointClassName)
            if klass:
                trigPointFinderClasses.append(klass)
            else:
                print >>sys.stderr, (
                    'Could not find trig point finder class %r! Has that '
                    'class been renamed or removed?' % trigPointClassName)
                sys.exit(1)

        database = Database(landmarkFinderClasses, trigPointFinderClasses,
                            limitPerLandmark=state['limitPerLandmark'],
                            maxDistance=state['maxDistance'],
                            minDistance=state['minDistance'])

        # Monkey-patch the new database instance to restore its state.
        for attr in ('_checksum', 'd', 'subjectCount', 'totalResidues',
                     'totalCoveredResidues', 'subjectInfo'):
            setattr(database, attr, state[attr])

        return database

    def checksum(self):
        """
        Compute a checksum for the specification and contents of this database.

        @return: A C{str} checksum.
        """
        # If the checksum is already known, return it.
        if self._checksum is not None:
            return self._checksum

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
        update(self.minDistance)
        update(self.bucketFactor)

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

        self._checksum = result.hexdigest()
        return self._checksum

    def checksumDirty(self):
        """
        Check if the database has an up-to-date checksum stored.

        @return: A C{bool}, C{True} if the in-memory checksum is valid,
            C{False} if not.
        """
        # This function exists for testing purposes.
        return self._checksum is None
