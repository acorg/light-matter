import logging
from io import StringIO
from collections import defaultdict, OrderedDict
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from Bio.File import as_handle

from light.subject import Subject, SubjectStore
from light.distance import scale
from light.checksum import Checksum
from light.parameters import Parameters
from light.reads import ScannedRead
from light.landmarks import landmarkNameFromHashkey
from light.trig import trigNameFromHashkey
from light.string import MultilineString


class Backend:
    """
    Maintain a collection of hashes for sequences (aka subjects) and provide
    for database index creation and search operations on them.

    @param subjectStore: A C{SubjectStore} instance, or C{None} if a
        C{SubjectStore} should be created.
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """

    # Implementation note: There are various methods in this class that
    # access self.params. But self.params is not put in place until the
    # 'configure' method has been called. We could add a 'checkConfigured'
    # method and call it on each function call if we wanted to be very
    # conservative and to supply good error messages if that happens.
    # Instead, we rely on the code that creates a Backend instance to
    # initialize it before trying to use it. If they do not, an
    # AttributeError will be raised on an attempted use of self.params
    # below. Hopefully this comment will indicate what's going wrong.

    DEFAULT_NAME = 'backend'
    SAVE_SUFFIX = '.lmbe'

    def __init__(self, subjectStore=None, filePrefix=None):
        self._configured = False
        self._subjectStore = subjectStore or SubjectStore()
        for method in ('getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects'):
            setattr(self, method, getattr(self._subjectStore, method))

        self._filePrefix = filePrefix
        self._totalCoveredResidues = 0

        # When you read the code in addSubject (below), it may look like
        # self.d could be a defaultdict(list). That does not work, though,
        # because a database JSON save followed by a load restores the
        # defaultdict as a normal dict.
        self.d = {}

    def configure(self, params, suggestedName=None, suggestedChecksum=None):
        """
        Configure this backend.

        @param params: A C{Parameters} instance.
        @param suggestedName: The C{str} suggested name for this backend. If
            the backend has already been configured (from a file restore) with
            a different name, the suggested name is ignored.
        @param suggestedChecksum: The C{int} suggested checksum for the
            backend, or C{None} if there is no initial value.
        @return: A 2-tuple consisting of the C{str} name of the backend and
            its checksum. These will either be the suggested values or those
            that were already in use (if the backend was already configured).
        """
        if self._configured:
            # We're already configured, possibly because we were created
            # via 'restore' or possibly because the WAMP connector has been
            # restarted. Check that the passed parameters match what we
            # already have in place.
            if params.compare(self.params) is not None:
                fp = StringIO()
                self.params._print(fp)
                original = fp.getvalue()
                fp = StringIO()
                params._print(fp)
                new = fp.getvalue()
                error = ('Already configured backend passed different '
                         'parameters.\nOriginal parameters\n%s\nLater '
                         'parameters\n%s' % (original, new))
                logging.warning(error)

            # There is no point checking the suggested name, as it's not an
            # error to be passed a suggested name that's different from the
            # one we're already using. This happens when a backend restores
            # itself from a file (which results in a call to configure) and
            # then reconnects to the WAMP server. The WAMP server will make
            # up a new suggested name and checksum because it cannot
            # recognize a previously connected backend just from the
            # ephemeral session id it receives when a connection is made.

            logging.debug('configure called again, with correct parameters')
        else:
            self.params = params
            self.name = suggestedName or self.DEFAULT_NAME

            if suggestedChecksum is None:
                self._checksum = self.initialChecksum(params, self.name)
            else:
                self._checksum = Checksum(suggestedChecksum)

            self._configured = True

        return self.name, self.checksum(), self.subjectCount()

    @staticmethod
    def initialChecksum(params, name):
        """
        Calculate an inital checksum for a backend.

        @param params: A C{Parameters} instance.
        @param name: The C{str} name of thd backend.
        @return: An initialized C{Checksum} instance for the given parameters
            and name.
        """
        return Checksum(params.checksum).update([name])

    def checksum(self):
        """
        Get the current checksum.

        @return: An C{int} checksum.
        """
        return self._checksum.value

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
            hash count for the subject.
        """
        subject = Subject(subject.id, subject.sequence, 0, subject.quality)
        preExisting, subjectIndex = self._subjectStore.add(subject,
                                                           subjectIndex)

        if preExisting:
            return (
                True, subjectIndex,
                self._subjectStore.getSubjectByIndex(subjectIndex).hashCount)

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

        return False, subjectIndex, subject.hashCount

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

    def find(self, read, storeFullAnalysis=False):
        """
        Given a read, compute all hashes for it, look up matching hashes and
        check which database sequences it matches.

        @param read: A C{dark.reads.AARead} instance.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A tuple of length three, containing
            1. Matches, a C{dict} keyed by subject index and whose values are
               as shown below.
            2. The number of hashes found in the read.
            3. A C{dict} of non-matching hashes, keyed by hash with values as
               returned by self.getHashes.
        """
        matches = defaultdict(list)
        nonMatchingHashes = {}
        readHashCount = 0
        scannedRead = self.scan(read)

        for readHash, hashInfo in self.getHashes(scannedRead).items():
            # Note that readHashCount is incremented for every hash in the
            # read, even those that are not found in the database. Basing
            # significance (in result.py) on a fraction of that overall
            # read hash count therefore takes into account the fact that
            # some hashes may have been missed. We may later want to do
            # that at a finer level of granularity, though.  E.g., by
            # considering where in the read the misses were.
            readHashCount += len(hashInfo['offsets'])
            try:
                subjectDict = self.d[readHash]
            except KeyError:
                # A hash that's in the read but not in our database.
                #
                # Note: this overwrites any existing value already stored
                #       in nonMatchingHashes[readHash]. Should the value be a
                #       list, or do we really not care? See
                #       https://github.com/acorg/light-matter/issues/301
                if storeFullAnalysis:
                    nonMatchingHashes[readHash] = hashInfo
            else:
                for (subjectIndex, subjectOffsets) in subjectDict.items():
                    matches[subjectIndex].append({
                        'landmark': hashInfo['landmark'],
                        'queryOffsets': hashInfo['offsets'],
                        'subjectOffsets': subjectOffsets,
                        'trigPoint': hashInfo['trigPoint']})

        return matches, readHashCount, nonMatchingHashes

    def shutdown(self, save, filePrefix):
        """
        Shut down the backend.

        @param save: If C{True}, save the backend state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        if save:
            self.save(filePrefix)

    def save(self, fpOrFilePrefix=None):
        """
        Save the backend database to a file.

        @param fpOrFilePrefix: A file pointer, or the C{str} prefix of a file
            name, or C{None}. If a C{str}, self.SAVE_SUFFIX is appended to get
            the full file name. If C{None}, self._filePrefix will be used as a
            file prefix unless it is also C{None}.
        @raises ValueError: If C{fpOrFilePrefix} and C{self._filePrefix} are
            both C{None}
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + self.SAVE_SUFFIX
        elif fpOrFilePrefix is None:
            if self._filePrefix is None:
                raise ValueError('save must be given an argument (or the '
                                 'database must have been restored from a '
                                 'file).')
            else:
                saveFile = self._filePrefix + self.SAVE_SUFFIX
        else:
            saveFile = fpOrFilePrefix

        state = {
            'name': self.name,
            'checksum': self.checksum(),
            'd': self.d,
            '_totalCoveredResidues': self._totalCoveredResidues,
        }

        with as_handle(saveFile, 'w') as fp:
            self.params.save(fp)
            self._subjectStore.save(fp)
            dump(state, fp)
            fp.write('\n')

    @classmethod
    def restore(cls, fpOrFilePrefix):
        """
        Load a database backend from a file.

        @param fpOrFilePrefix: A file pointer or the C{str} prefix of a file
            name. If a C{str}, self.SAVE_SUFFIX is appended to get the full
            file name.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved database backend file. Or if valid JSON
            cannot be loaded from C{fp}.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + cls.SAVE_SUFFIX
            filePrefix = fpOrFilePrefix
        else:
            saveFile = fpOrFilePrefix
            filePrefix = None

        with as_handle(saveFile) as fp:
            params = Parameters.restore(fp)
            subjectStore = SubjectStore.restore(fp)
            state = loads(fp.readline()[:-1])

        new = cls(subjectStore, filePrefix=filePrefix)
        new.configure(params, suggestedName=state['name'],
                      suggestedChecksum=state['checksum'])

        for attr in ('d', '_totalCoveredResidues'):
            setattr(new, attr, state[attr])

        return new

    def print_(self, margin=''):
        """
        Print information about the database index.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the parameters.
        """
        result = MultilineString(margin=margin)
        append = result.append
        indent = result.indent
        outdent = result.outdent

        result.extend([
            'Name: %s' % self.name,
            'Hash count: %d' % len(self.d),
            'Checksum: %s' % self.checksum(),
        ])

        append('Subjects (with offsets) by hash:')
        indent()
        landmarkCount = defaultdict(int)
        trigCount = defaultdict(int)
        for hash_, subjects in sorted(self.d.items()):
            append(hash_)
            # The split on ':' corresponds to the use of ':' above in
            # self.hash() to make a hash key.
            landmarkHashkey, trigHashkey, distance = hash_.split(':')
            landmarkCount[landmarkHashkey] += 1
            trigCount[trigHashkey] += 1
            indent()
            for subjectIndex, offsets in sorted(subjects.items()):
                append('%s %r' % (subjectIndex, offsets))
            outdent()

        outdent()
        append('Landmark symbol counts:')
        indent()

        for hash_, count in sorted(landmarkCount.items()):
            append('%s (%s): %d' % (
                landmarkNameFromHashkey(hash_), hash_, count))

        outdent()
        append('Trig point symbol counts:')
        indent()

        for hash_, count in sorted(trigCount.items()):
            append('%s (%s): %d' % (
                trigNameFromHashkey(hash_), hash_, count))

        return str(result)
