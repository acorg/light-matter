import sys
from collections import defaultdict, OrderedDict
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from light.subject import Subject, SubjectStore
from light.distance import scale
from light.checksum import Checksum
from light.parameters import Parameters
from light.reads import ScannedRead
from light.landmarks import landmarkNameFromHashkey
from light.trig import trigNameFromHashkey


class Backend:
    """
    Maintain a collection of hashes for sequences (aka subjects) and provide
    for database index creation and search operations on them.

    @param params: A C{Parameters} instance.
    @param name: The C{str} name of this backend.
    @param checksum: The C{int} initial checksum for the backend, or C{None}
        if there is no initial value.
    @param subjectStore: A C{SubjectStore} instance, or C{None} if a
        C{SubjectStore} should be created.
    """
    DEFAULT_NAME = 'backend'

    def __init__(self, params, name=DEFAULT_NAME, checksum=None,
                 subjectStore=None):
        self.params = params
        self.name = name
        self._subjectStore = subjectStore or SubjectStore()
        for method in ('getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects'):
            setattr(self, method, getattr(self._subjectStore, method))

        # If no initial checksum value was given, initialize the checksum
        # with the checksum from our parameters and our name.
        if checksum is None:
            self._checksum = Checksum(params.checksum).update([name])
        else:
            self._checksum = Checksum(checksum)

        self._totalCoveredResidues = 0

        # When you read the code in addSubject (below), it may look like
        # self.d could be a defaultdict(list). That does not work, though,
        # because a database JSON save followed by a load restores the
        # defaultdict as a normal dict.
        self.d = {}

    def checksum(self):
        """
        Get the current checksum.

        @return: An C{int} checksum.
        """
        return self._checksum.checksum

    def subjectCount(self):
        """
        How many subjects are stored in this backend?

        @return: An C{int} number of subjects.
        """
        return len(self._subjectStore)

    def hashCount(self):
        """
        How many hashes are in the backend database?

        @return: An C{int} number of hashes.
        """
        return len(self.d)

    def totalResidues(self):
        """
        How many AA residues are there in the subjects stored in this backend?

        @return: An C{int} number of residues.
        """
        return sum(len(s) for s in self._subjectStore.getSubjects())

    def totalCoveredResidues(self):
        """
        How many AA residues are covered by landmarks and trig points in the
        subjects stored in this backend?

        @return: The C{int} number of residues that are covered.
        """
        return self._totalCoveredResidues

    def scan(self, sequence):
        """
        Make an instance of C{light.reads.ScannedRead} from a sequence.

        @param sequence: a C{dark.read.AARead} instance.
        @return: a C{light.reads.ScannedRead} instance.
        """
        scannedSequence = ScannedRead(sequence)

        append = scannedSequence.landmarks.append
        for landmarkFinder in self.params.landmarkFinders:
            for landmark in landmarkFinder.find(sequence):
                append(landmark)

        append = scannedSequence.trigPoints.append
        for trigFinder in self.params.trigPointFinders:
            for trigPoint in trigFinder.find(sequence):
                append(trigPoint)

        return scannedSequence

    def hash(self, landmark, trigPoint):
        """
        Compute a hash key to store information about a landmark / trig point
        association for a read.

        @param landmark: A C{light.features.Landmark} instance.
        @param trigPoint: A C{light.features.TrigPoint} instance.
        @return: A C{str} hash key based on the landmark, the trig point,
            and the distance between them.
        """
        distance = scale(trigPoint.offset - landmark.offset,
                         self.params.distanceBase)
        return '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                             distance)

    def getScannedPairs(self, scannedSequence):
        """
        Get the (landmark, trigPoint) pairs from a ScannedRead instance.

        @param scannedSequence: A C{light.reads.ScannedRead} instance.
        @return: A generator yielding (landmark, trigPoint) pairs, as returned
            by C{light.reads.ScannedRead.getPairs}.
        """
        return scannedSequence.getPairs(
            limitPerLandmark=self.params.limitPerLandmark,
            maxDistance=self.params.maxDistance,
            minDistance=self.params.minDistance)

    def getHashes(self, scannedSequence):
        """
        Get all (landmark, trigPoint) hashes from a scanned sequence and
        collect the offsets at which the (landmark, trigPoint) pair occurs.

        @param scannedSequence: A C{light.reads.ScannedRead} instance.
        @return: A C{dict} keyed by (landmark, trigPoint) hash, whose values
            are C{dict}s containing the first (landmark, trigPoint) pair that
            generated that hash, and a list of all offsets into the read
            where the (landmark, trigPoint) pair was found.
        """
        # For testing reasons, use an ordered dict to hold hash information.
        # Our database and result code do not need the dict to be ordered.
        # But a deterministic order of hashes makes it simple to write reliable
        # tests. If we ever get really serious about speed we may want to use
        # a regular dict instead and make the tests do more digging / sorting
        # in results.
        hashes = OrderedDict()

        for (landmark, trigPoint) in self.getScannedPairs(scannedSequence):
            hash_ = self.hash(landmark, trigPoint)
            try:
                hashInfo = hashes[hash_]
            except KeyError:
                hashes[hash_] = {
                    'landmark': landmark,
                    'offsets': [[landmark.offset, trigPoint.offset]],
                    'trigPoint': trigPoint,
                }
            else:
                hashInfo['offsets'].append([landmark.offset, trigPoint.offset])

        return hashes

    def addSubject(self, subject, subjectIndex=None):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: A C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @param subjectIndex: A C{str} representing the index of the subject as
            known by the database front end. If C{None} the SubjectStore we
            call addSubject on will assign an index.
        @return: A tuple of 1) a C{bool} to indicate whether the subject was
            already in the database, 2) the C{str} subject index, and 3) the
            C{str} name of this backend.
        """
        subject = Subject(subject.id, subject.sequence, 0, subject.quality)
        preExisting, subjectIndex = self._subjectStore.add(subject,
                                                           subjectIndex)

        if preExisting:
            return True, subjectIndex, self.name

        scannedSubject = self.scan(subject)
        self._totalCoveredResidues += len(scannedSubject.coveredIndices())

        for landmark, trigPoint in self.getScannedPairs(scannedSubject):
            hash_ = self.hash(landmark, trigPoint)
            subject.hashCount += 1
            try:
                subjectDict = self.d[hash_]
            except KeyError:
                self.d[hash_] = subjectDict = {}

            # Don't use a tuple for the offsets because JSON save/load will
            # convert it to a list and we'll need to convert all those lists
            # to tuples on database load.
            offsets = [landmark.offset, trigPoint.offset]
            try:
                subjectDict[subjectIndex].append(offsets)
            except KeyError:
                subjectDict[subjectIndex] = [offsets]

        self._checksum.update((subject.id, subject.sequence))

        return False, subjectIndex, self.name

    def find(self, read, significanceMethod, scoreMethod,
             significanceFraction, storeFullAnalysis):
        """
        Given a read, compute all hashes for it, look up matching hashes and
        check which database sequences it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The C{str} name of the method used to
            calculate which histogram bins are considered significant.
        @param scoreMethod: The C{str} name of the method used to calculate the
            score of a bin which is considered significant.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A list of length three, containing
            1. Matches, a C{dict} keyed by subject index and whose values are
               as shown below.
            2. The number of hashes found in the read.
            3. A C{dict} of non-matching hashes, keyed by hash with values as
               returned by self.getHashes.
        """
        matches = defaultdict(list)
        nonMatchingHashes = {}
        hashCount = 0
        scannedRead = self.scan(read)

        for hash_, hashInfo in self.getHashes(scannedRead).items():
            # Note that hashCount is incremented for every hash, even those
            # that are not in the database. Basing significance (in
            # result.py) on a fraction of that overall count therefore
            # takes into account the fact that some hashes may have been
            # missed. We may want to do that at a finer level of
            # granularity, though.  E.g., by considering where in the read
            # the misses were.
            hashCount += len(hashInfo['offsets'])
            try:
                subjectDict = self.d[hash_]
            except KeyError:
                # A hash that's in the read but not in our database.
                if storeFullAnalysis:
                    nonMatchingHashes[hash_] = hashInfo
            else:
                for (subjectIndex, subjectOffsets) in subjectDict.items():
                    matches[subjectIndex].append({
                        'landmark': hashInfo['landmark'],
                        'queryOffsets': hashInfo['offsets'],
                        'subjectOffsets': subjectOffsets,
                        'trigPoint': hashInfo['trigPoint']})

        return matches, hashCount, nonMatchingHashes

    def save(self, fp=sys.stdout):
        """
        Save the backend database to a file.

        @param fp: A file pointer.
        """
        self.params.save(fp)
        self._subjectStore.save(fp)
        state = {
            'name': self.name,
            'checksum': self.checksum(),
            'd': self.d,
            '_totalCoveredResidues': self._totalCoveredResidues,
        }
        dump(state, fp)
        fp.write('\n')

    @classmethod
    def restore(cls, fp=sys.stdin):
        """
        Load a database backend from a file.

        @param fp: A file pointer.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved database backend file. Or if valid JSON
            cannot be loaded from C{fp}.
        """
        params = Parameters.restore(fp)
        subjectStore = SubjectStore.restore(fp)

        state = loads(fp.readline()[:-1])
        new = cls(params, name=state['name'], checksum=state['checksum'],
                  subjectStore=subjectStore)

        for attr in ('d', '_totalCoveredResidues'):
            setattr(new, attr, state[attr])

        return new

    def print_(self, fp=sys.stdout):
        """
        Print information about the database index.

        @param fp: A file pointer to write to.
        """
        print('Name: %s' % self.name, file=fp)
        print('Hash count: %d' % len(self.d), file=fp)
        print('Checksum: %s' % self.checksum(), file=fp)

        print('Subjects (with offsets) by hash:', file=fp)
        landmarkCount = defaultdict(int)
        trigCount = defaultdict(int)
        for hash_, subjects in sorted(self.d.items()):
            print('  ', hash_, file=fp)
            # The split on ':' corresponds to the use of ':' above in
            # self.hash() to make a hash key.
            landmarkHashkey, trigHashkey, distance = hash_.split(':')
            landmarkCount[landmarkHashkey] += 1
            trigCount[trigHashkey] += 1
            for subjectIndex, offsets in sorted(subjects.items()):
                print('    %s %r' % (subjectIndex, offsets), file=fp)

        print('Landmark symbol counts:', file=fp)
        for hash_, count in sorted(landmarkCount.items()):
            print('  %s (%s): %d' % (
                landmarkNameFromHashkey(hash_), hash_, count), file=fp)

        print('Trig point symbol counts:', file=fp)
        for hash_, count in sorted(trigCount.items()):
            print('  %s (%s): %d' % (
                trigNameFromHashkey(hash_), hash_, count), file=fp)
