import sys
from collections import defaultdict
from ujson import dump, dumps, load
from binascii import crc32
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
    SIGNIFICANCE_FRACTION_DEFAULT = 0.25

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
        if bucketFactor <= 0:
            raise ValueError('bucketFactor must be > 0.')
        else:
            self.bucketFactor = bucketFactor
        # It may look like self.d should be a defaultdict(list). But that
        # will not work because a database JSON save followed by a load
        # will restore the defaultdict as a vanilla dict.
        self.d = {}
        self.subjectCount = 0
        self.totalResidues = 0
        self.totalCoveredResidues = 0
        self.subjectInfo = []
        # Create instances of the landmark and trig point finder classes.
        self.landmarkFinders = []
        for landmarkFinderClass in self.landmarkFinderClasses:
            self.landmarkFinders.append(landmarkFinderClass())
        self.trigPointFinders = []
        for trigPointFinderClass in self.trigPointFinderClasses:
            self.trigPointFinders.append(trigPointFinderClass())
        self._initializeChecksum()

    def _initializeChecksum(self):
        """
        Set the initial checksum, based on the database parameters.
        """
        self.checksum = 0x0  # An arbitrary starting checksum.
        key = attrgetter('NAME')
        self._updateChecksum(
            [f.NAME for f in sorted(self.landmarkFinders, key=key)] +
            [f.NAME for f in sorted(self.trigPointFinders, key=key)] +
            map(str, (self.limitPerLandmark, self.maxDistance,
                      self.minDistance, self.bucketFactor)))

    def _updateChecksum(self, strings):
        """
        Update the checksum for this database.

        @param strings: A C{list} of strings to update the current checksum
            with.
        """
        update = '\0'.join(strings) + '\0'
        self.checksum = crc32(update, self.checksum) & 0xFFFFFFFF

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
        subjectInfo = (subject.id, subject.sequence)
        self._updateChecksum(subjectInfo)
        self.subjectInfo.append(subjectInfo)
        scannedSubject = ScannedSubject(subject)
        subjectIndex = self.subjectCount
        self.subjectCount += 1
        self.totalResidues += len(subject)

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
            try:
                subjectDict = self.d[key]
            except KeyError:
                self.d[key] = subjectDict = {}

            try:
                subjectDict[str(subjectIndex)].append(landmark.offset)
            except KeyError:
                subjectDict[str(subjectIndex)] = [landmark.offset]

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
            'checksum': self.checksum,
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

    def find(self, read, significanceFraction=None, storeAnalysis=False):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if significanceFraction is None:
            significanceFraction = self.SIGNIFICANCE_FRACTION_DEFAULT

        scannedRead = ScannedRead(read)

        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(read):
                scannedRead.trigPoints.append(trigPoint)

        matches = defaultdict(list)
        hashCount = 0

        for landmark, trigPoint in scannedRead.getPairs(
                limitPerLandmark=self.limitPerLandmark,
                maxDistance=self.maxDistance, minDistance=self.minDistance):
            hashCount += 1
            key = self.key(landmark, trigPoint)
            try:
                subjectDict = self.d[key]
            except KeyError:
                # A hash that's in the read but not in our database. We
                # should eventually keep these mismatches and use them to
                # help determine how good a match against a subject is.
                #
                # Note that hashCount is incremented for every pair, even
                # ones that are not in the database. Basing significance on
                # a fraction of that overall count does take into account
                # the fact that some hashes may have been missed. We may
                # want to do that at a finer level of granularity, though.
                # E.g., consider _where_ in the read the misses were.
                pass
            else:
                for subjectIndex, subjectOffsets in subjectDict.iteritems():
                    subjectIndex = int(subjectIndex)
                    subjectLength = len(self.subjectInfo[subjectIndex][1])
                    matches[subjectIndex].append({
                        'landmarkLength': landmark.length,
                        'landmarkName': landmark.name,
                        'readOffset': landmark.offset,
                        'subjectLength': subjectLength,
                        'subjectOffsets': subjectOffsets,
                        'trigPointName': trigPoint.name,
                    })

        return Result(read, matches, hashCount, significanceFraction,
                      self.bucketFactor, storeAnalysis=storeAnalysis)

    def save(self, fp=sys.stdout):
        """
        Save the database to a file.

        @param fp: A file pointer.
        """
        state = {
            'checksum': self.checksum,
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
                            minDistance=state['minDistance'],
                            bucketFactor=state['bucketFactor'])

        # Monkey-patch the new database instance to restore its state.
        for attr in ('checksum', 'd', 'subjectCount', 'totalResidues',
                     'totalCoveredResidues', 'subjectInfo'):
            setattr(database, attr, state[attr])

        return database
