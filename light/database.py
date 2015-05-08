import argparse
import sys
from warnings import warn
from collections import defaultdict, OrderedDict
try:
    from ujson import dump, load
except ImportError:
    from json import dump, load

from dark.fasta import combineReads, FastaReads
from dark.reads import AARead

from light.distance import scale
from light.checksum import Checksum
from light.subject import Subject
from light.parameters import Parameters
from light.reads import ScannedRead
from light.result import Result
from light.landmarks import (
    findLandmarks, landmarkNameFromHashkey, ALL_LANDMARK_CLASSES,
    DEFAULT_LANDMARK_CLASSES)
from light.trig import (
    findTrigPoints, trigNameFromHashkey, ALL_TRIG_CLASSES,
    DEFAULT_TRIG_CLASSES)


class _DatabaseMixin(object):
    """
    Methods and data that can be used by databases and their back-ends.

    @param params: A C{Parameters} instance.
    """
    def __init__(self, params):
        self.params = params
        self.subjectCount = 0
        # Base the database checksum on the parameters checksum.
        self._checksum = Checksum(params.checksum)

    def _getChecksum(self):
        """
        Get the current checksum.

        @return: An C{int} checksum.
        """
        return self._checksum.checksum

    def _setChecksum(self, value):
        """
        Set the current checksum.

        @param value: An C{int} value to set the checksum to.
        """
        self._checksum.checksum = value

    checksum = property(_getChecksum, _setChecksum)

    def scan(self, sequence):
        """
        Make an instance of C{light.reads.ScannedRead} from a sequence.

        @param sequence: a C{dark.read.AARead} instance.
        @return: a C{light.reads.ScannedRead} instance.
        """
        scannedSequence = ScannedRead(sequence)
        for landmarkFinder in self.params.landmarkFinders:
            for landmark in landmarkFinder.find(sequence):
                scannedSequence.landmarks.append(landmark)

        for trigFinder in self.params.trigPointFinders:
            for trigPoint in trigFinder.find(sequence):
                scannedSequence.trigPoints.append(trigPoint)
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
                    'landmarkOffsets': [landmark.offset],
                    'trigPointOffsets': [trigPoint.offset],
                    'trigPoint': trigPoint,
                }
            else:
                hashInfo['trigPointOffsets'].append(trigPoint.offset)
                hashInfo['landmarkOffsets'].append(landmark.offset)

        return hashes


class Database(_DatabaseMixin):
    """
    Maintain a collection of sequences ("subjects") and provide for database
    insertion and look-up operations on them.

    @param params: A C{Parameters} instance.
    """
    def __init__(self, params, backendConnector=None):
        super().__init__(params)
        self.totalResidues = 0
        self.totalCoveredResidues = 0
        self.hashCount = 0
        self._subjectInfo = []
        self._idSequenceCache = {}
        self._backendConnector = backendConnector or SimpleConnector(
            Backend(params))
        self._backendSubjectIndex = {}

    def __str__(self):
        return '%s: %d sequences, %d residues, %.2f%% coverage' % (
            self.__class__.__name__, self.subjectCount, self.totalResidues,
            float(self.totalCoveredResidues) / self.totalResidues * 100.0)

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @return: The C{int} subject index of the added subject.
        """
        try:
            return self._idSequenceCache[subject]
        except KeyError:
            pass

        subjectIndex = self.subjectCount
        self.subjectCount += 1
        self._idSequenceCache[subject] = subjectIndex
        self.totalResidues += len(subject)

        backendName, backendSubjectIndex, coveredResidues, hashCount = \
            self._backendConnector.addSubject(subject)

        self.hashCount += hashCount

        try:
            indexMap = self._backendSubjectIndex[backendName]
        except KeyError:
            indexMap = self._backendSubjectIndex[backendName] = {}

        indexMap[backendSubjectIndex] = subjectIndex

        self.totalCoveredResidues += coveredResidues
        self._checksum.update((subject.id, subject.sequence, str(hashCount)))
        self._subjectInfo.append((subject.id, subject.sequence, hashCount))

        return subjectIndex

    def getSubject(self, subject):
        """
        Return information about a subject, given its index in the database.
        Or, do the reverse, given a subject return its index.

        @param subject: Either an C{int} subject index or a C{dark.read.AARead}
            representing a subject whose index is to be looked up.
        @return: If an index is passed, return a C{Subject} instance. If an
            C{AARead} subject is passed, return its index.
        @raises IndexError: if an C{int} index is passed and it is invalid.
        @raises KeyError: if an C{AARead} subject is passed and it is not in
            the database.
        """
        if isinstance(subject, int):
            return Subject(*self._subjectInfo[subject])
        else:
            try:
                return self._idSequenceCache[subject]
            except KeyError:
                # Be user friendly and raise a key error containing the
                # subject id, instead of the md5 sum of the id and
                # sequence.
                raise KeyError(subject.id)

    def getSubjects(self):
        """
        Return information about all database subjects.

        @return: a generator that yields C{Subject} instances.
        """
        return (self.getSubject(i) for i in range(self.subjectCount))

    def find(self, read, significanceMethod=None, scoreMethod=None,
             significanceFraction=None, storeFullAnalysis=False):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The name of the method used to calculate
            which histogram bins are considered significant.
        @param scoreMethod: The C{str} name of the method used to calculate the
            score of a bin which is considered significant.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if significanceMethod is None:
            significanceMethod = self.params.DEFAULT_SIGNIFICANCE_METHOD
        if scoreMethod is None:
            scoreMethod = self.params.DEFAULT_SCORE_METHOD
        if significanceFraction is None:
            significanceFraction = self.params.DEFAULT_SIGNIFICANCE_FRACTION

        allMatches = defaultdict(list)
        allNonMatchingHashes = {}
        hashCount = 0

        for result in self._backendConnector.find(
                read, significanceMethod, scoreMethod,
                significanceFraction, storeFullAnalysis):
            backendName, matches, hashCount, nonMatchingHashes = result
            # TODO: This is probably wrong...
            for nonMatchingHash in nonMatchingHashes:
                if nonMatchingHash not in allNonMatchingHashes:
                    allNonMatchingHashes[nonMatchingHash] = nonMatchingHashes[
                        nonMatchingHash]
            for backendSubjectIndex in matches:
                # Make sure we have not have seen this subject before. If
                # we have, it would mean that two backends are reporting
                # results for the same subject. We may allow that later,
                # but for now it should be an error.
                assert(backendSubjectIndex not in allMatches)

                # Convert the backend subject index to one of our indices.
                subjectIndex = self._backendSubjectIndex[backendName][
                    backendSubjectIndex]
                allMatches[subjectIndex] = matches[backendSubjectIndex]

        return Result(self.scan(read), allMatches, hashCount,
                      significanceMethod, scoreMethod, significanceFraction,
                      self, allNonMatchingHashes,
                      storeFullAnalysis=storeFullAnalysis)

    def save(self, fp=sys.stdout):
        """
        Save the database state to a file.

        @param fp: A file pointer.
        """
        self.params.save(fp)
        state = {
            'checksum': self.checksum,
            'subjectCount': self.subjectCount,
            'totalResidues': self.totalResidues,
            'totalCoveredResidues': self.totalCoveredResidues,
            '_subjectInfo': self._subjectInfo,
            '_backendConnectorClass':
                self._backendConnector.__class__.__name__,
            '_backendSubjectIndex': self._backendSubjectIndex,
        }
        dump(state, fp)
        fp.write('\n')

    @staticmethod
    def restore(fp=sys.stdin):
        """
        Load a database from a file.

        @param fp: A file pointer.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved database file. Or if valid JSON cannot be
            loaded from C{fp}.
        """
        params = Parameters.restore(fp)
        backend = Backend(params)
        state = load(fp)

        if state['_backendConnectorClass'] == 'SimpleConnector':
            backendConnector = SimpleConnector(backend)
        else:
            raise ValueError('Unknown backend connector class %r.' %
                             state['_backendConnectorClass'])
        database = Database(params, backendConnector)

        # Monkey-patch the new database instance to restore the rest of its
        # state.
        for attr in ('subjectCount', 'totalResidues', 'totalCoveredResidues',
                     '_backendSubjectIndex'):
            setattr(database, attr, state[attr])

        # Re-initialize the database checksum.
        database.checksum = state['checksum']

        # The _subjectInfo entries were originally tuples, but JSON I/O turns
        # these into lists. For safety, convert them back to tuples.
        _subjectInfo = state['_subjectInfo']
        for index, subjectInfo in enumerate(_subjectInfo):
            _subjectInfo[index] = tuple(subjectInfo)
        database._subjectInfo = _subjectInfo

        # Restore the id/sequence cache.
        getSubject = database.getSubject
        for subjectIndex in range(database.subjectCount):
            database._idSequenceCache[getSubject(subjectIndex)] = subjectIndex

        return database

    def print_(self, fp=sys.stdout, printHashes=False):
        """
        Print information about the database.

        @param fp: A file pointer to write to.
        @param printHashes: If C{True}, print all hashes and associated
            subjects.
        """
        self.params.print_(fp)
        try:
            coverage = (float(self.totalCoveredResidues) /
                        self.totalResidues * 100.0)
        except ZeroDivisionError:
            coverage = 0.0

        print('Subject count: %s' % self.subjectCount, file=fp)
        print('Hash count: %s' % self.hashCount, file=fp)
        print('Total residues: %d' % self.totalResidues, file=fp)
        print('Coverage: %.2f%%' % coverage, file=fp)
        print('Checksum: %s' % self.checksum, file=fp)

        if printHashes and self.subjectCount:
            self._backendConnector.print_(fp)

    def emptyCopy(self):
        """
        Copy the current database, with no subjects.

        @return: A new L{light.database.Database} instance that no subjects
            have been added to.
        """
        return self.__class__(self.params, self._backendConnector.emptyCopy())


class Backend(_DatabaseMixin):
    """
    Maintain a collection of hashes for sequences ("subjects") and provide for
    database index creation and search operations on them.

    @param params: A C{Parameters} instance.
    """
    def __init__(self, params):
        super().__init__(params)

        # It may look like self.d should be a defaultdict(list). But that
        # will not work because a database JSON save followed by a load
        # will restore the defaultdict as a vanilla dict.
        self.d = {}

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @return: A list of length three, containing
            1. The C{str} subject index of the added subject.
            2. The number of residues covered in the subject by its landmarks
               and trig points.
            3. The number of hashes in the subject (a hash is a landmark /
               trigpoint / distance combination).
        """
        subjectIndex = str(self.subjectCount)
        self.subjectCount += 1
        scannedSubject = self.scan(subject)
        coveredResidues = len(scannedSubject.coveredIndices())
        hashCount = 0

        for landmark, trigPoint in self.getScannedPairs(scannedSubject):
            hash_ = self.hash(landmark, trigPoint)
            hashCount += 1
            try:
                subjectDict = self.d[hash_]
            except KeyError:
                self.d[hash_] = subjectDict = {}

            try:
                subjectDict[subjectIndex].append(landmark.offset)
            except KeyError:
                subjectDict[subjectIndex] = [landmark.offset]

        self._checksum.update((subject.id, subject.sequence, str(hashCount)))

        return subjectIndex, coveredResidues, hashCount

    def find(self, read, significanceMethod, scoreMethod, significanceFraction,
             storeFullAnalysis):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The name of the method used to calculate
            which histogram bins are considered significant.
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
            hashCount += len(hashInfo['landmarkOffsets'])
            try:
                subjectDict = self.d[hash_]
            except KeyError:
                # A hash that's in the read but not in our database.
                if storeFullAnalysis:
                    nonMatchingHashes[hash_] = hashInfo
            else:
                for (subjectIndex,
                     subjectLandmarkOffsets) in subjectDict.items():
                    matches[subjectIndex].append({
                        'landmark': hashInfo['landmark'],
                        'queryLandmarkOffsets': hashInfo['landmarkOffsets'],
                        'queryTrigPointOffsets': hashInfo['trigPointOffsets'],
                        'subjectLandmarkOffsets': subjectLandmarkOffsets,
                        'trigPoint': hashInfo['trigPoint'],
                    })

        return matches, hashCount, nonMatchingHashes

    def save(self, fp=sys.stdout):
        """
        Save the database to a file.

        @param fp: A file pointer.
        """
        self.params.save(fp)
        state = {
            'checksum': self.checksum,
            'd': self.d,
            'subjectCount': self.subjectCount,
        }
        dump(state, fp)

    @staticmethod
    def restore(fp=sys.stdin):
        """
        Load a database backend from a file.

        @param fp: A file pointer.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved database backend file. Or if valid JSON
            cannot be loaded from C{fp}.
        """
        params = Parameters.restore(fp)
        backend = Backend(params)
        state = load(fp)

        for attr in ('checksum', 'd', 'subjectCount'):
            setattr(backend, attr, state[attr])

        return backend

    def print_(self, fp=sys.stdout):
        """
        Print information about the database index.

        @param fp: A file pointer to write to.
        """
        print('Hash count: %d' % len(self.d), file=fp)
        print('Checksum: %s' % self.checksum, file=fp)

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

    def emptyCopy(self):
        """
        Copy the current backend, with an empty index.

        @return: A new L{light.database.Backend} instance with an empty
            index.
        """
        return self.__class__(self.params)


class SimpleConnector():
    """
    Connect a Database instance to a single Backend instances.

    @param params: A C{Parameters} instance.
    """
    NAME = 'localhost'

    def __init__(self, backend):
        self._backend = backend

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @return: A tuple of length four, containing
            1. The name of this backend.
            2. The C{int} subject index of the added subject.
            3. The number of residues covered in the subject by its landmarks
               and trig points, and
            4. the number of hashes in the subject (a hash is a landmark /
               trigpoint / distance combination).
        """
        subjectIndex, coveredResidues, hashCount = self._backend.addSubject(
            subject)
        return (self.NAME, str(subjectIndex), coveredResidues, hashCount)

    def find(self, read, significanceMethod, scoreMethod, significanceFraction,
             storeFullAnalysis):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The name of the method used to calculate
            which histogram bins are considered significant.
        @param scoreMethod: The C{str} name of the method used to calculate the
            score of a bin which is considered significant.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{tuple} containing a single C{tuple} of length four,
            containing:
            1. The C{str} name of this backend.
            2. Matches, a C{dict} keyed by subject index and whose values are
               as shown below.
            3. The number of hashes found in the read.
            4. A C{dict} of non-matching hashes, keyed by hash with values as
               returned by self.getHashes.
        """
        matches, hashCount, nonMatchingHashes = self._backend.find(
            read, significanceMethod, scoreMethod, significanceFraction,
            storeFullAnalysis)
        return ((self.NAME, matches, hashCount, nonMatchingHashes),)

    def print_(self, fp=sys.stdout):
        """
        Print information about the database backends.

        @param fp: A file pointer to write to.
        """
        print('Backend 0:', file=fp)
        return self._backend.print_(fp)

    def emptyCopy(self):
        """
        Copy the current connector, with an empty index.

        @return: A new L{light.database.SimpleConnector} instance with a
            backend with an empty index.
        """
        return self.__class__(self._backend.emptyCopy())


class DatabaseSpecifier(object):
    """
    Helper class for creating or loading a database and populating it.

    @param allowCreation: If True, add options that permit the creation
        of a new database.
    @param allowPopulation: If True, add options that allow the automatic
        addition of sequences to the database.
    @param allowInMemory: If True, add the option that allows the user to
        pass an in-memory database (e.g., as constructed in ipython or
        iPythonNotebook or programmatically).
    @param allowFromFile: If True, add the option that allows the user to
        specify a pre-existing database.
    @raise ValueError: If the allow options do not permit creation or loading.
    """
    def __init__(self, allowCreation=True, allowPopulation=True,
                 allowInMemory=True, allowFromFile=True):
        if not (allowCreation or allowFromFile or allowInMemory):
            raise ValueError('You must either allow database creation, '
                             'loading a database from a file, or passing an '
                             'in-memory database.')
        self._allowCreation = allowCreation
        self._allowPopulation = allowPopulation
        self._allowInMemory = allowInMemory
        self._allowFromFile = allowFromFile

    def addArgsToParser(self, parser):
        """
        Add standard database creation or loading arguments to an argparse
        parser, depending on the allowable ways of specifying a database.

        @param parser: An C{argparse.ArgumentParser} instance.
        """
        if self._allowFromFile:
            parser.add_argument(
                '--databaseFile',
                help='The name of a file containing a saved database')

        if self._allowCreation:
            parser.add_argument(
                '--landmark', action='append', dest='landmarkFinderNames',
                choices=sorted(cl.NAME for cl in ALL_LANDMARK_CLASSES),
                help=('The name of a landmark finder to use. May be specified '
                      'multiple times.'))

            parser.add_argument(
                '--trig', action='append', dest='trigFinderNames',
                choices=sorted(cl.NAME for cl in ALL_TRIG_CLASSES),
                help=('The name of a trig point finder to use. May be '
                      'specified multiple times.'))

            parser.add_argument(
                '--defaultLandmarks', action='store_true', default=False,
                help=('If specified, use the default landmark finders: %s' %
                      sorted(cl.NAME for cl in
                             DEFAULT_LANDMARK_CLASSES)))

            parser.add_argument(
                '--defaultTrigPoints', action='store_true', default=False,
                help=('If specified, use the default trig point finders: %s' %
                      sorted(cl.NAME for cl in DEFAULT_TRIG_CLASSES)))

            parser.add_argument(
                '--limitPerLandmark', type=int,
                default=Parameters.DEFAULT_LIMIT_PER_LANDMARK,
                help=('A limit on the number of pairs to yield per landmark '
                      'per read.'))

            parser.add_argument(
                '--maxDistance', type=int,
                default=Parameters.DEFAULT_MAX_DISTANCE,
                help='The maximum distance permitted between yielded pairs.')

            parser.add_argument(
                '--minDistance', type=int,
                default=Parameters.DEFAULT_MIN_DISTANCE,
                help='The minimum distance permitted between yielded pairs.')

            parser.add_argument(
                '--distanceBase', type=float,
                default=Parameters.DEFAULT_DISTANCE_BASE,
                help=('The distance between a landmark and a trig point is '
                      'scaled to be its logarithm using this base. This '
                      'reduces sensitivity to relatively small differences in '
                      'distance.'))

        if self._allowPopulation:
            parser.add_argument(
                '--databaseFasta',
                help=('The name of a FASTA file containing the sequences that '
                      'should be added to the database.'))

            parser.add_argument(
                '--databaseSequence', action='append',
                dest='databaseSequences', metavar='"id sequence"',
                help=('Amino acid sequences to add to the database. The '
                      'sequence id will be the text up to the last space, if '
                      'any, otherwise will be automatically assigned.'))

    def getDatabaseFromArgs(self, args):
        """
        Read an existing database (if args.database is given) or create
        one from args.

        There is an order of preference in examining the arguments used to
        specify a database: pre-existing in a file (via --databaseFile),
        and then via the creation of a new database (many args). There are
        currently no checks to make sure the user isn't trying to do more
        than one of these at the same time (e.g., using both --databaseFile
        and --defaultLandmarks), the one with the highest priority is silently
        acted on first.

        @param args: Command line arguments as returned by the C{argparse}
            C{parse_args} method.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        if self._allowFromFile and args.databaseFile:
            with open(args.databaseFile) as fp:
                database = Database.restore(fp)

        elif self._allowCreation:
            landmarkClasses = (
                DEFAULT_LANDMARK_CLASSES if args.defaultLandmarks
                else findLandmarks(args.landmarkFinderNames))

            trigClasses = (
                DEFAULT_TRIG_CLASSES if args.defaultTrigPoints
                else findTrigPoints(args.trigFinderNames))

            if len(landmarkClasses) + len(trigClasses) == 0:
                warn("Creating a database with no landmark or trig point "
                     "finders. Hope you know what you're doing!")

            params = Parameters(landmarkClasses, trigClasses,
                                limitPerLandmark=args.limitPerLandmark,
                                maxDistance=args.maxDistance,
                                minDistance=args.minDistance,
                                distanceBase=args.distanceBase)
            backend = Backend(params)
            backendConnector = SimpleConnector(backend)
            database = Database(params, backendConnector)

        if self._allowPopulation:
            for read in combineReads(args.databaseFasta,
                                     args.databaseSequences, readClass=AARead):
                database.addSubject(read)

        return database

    def getDatabaseFromKeywords(
            self, landmarkNames=None, trigPointNames=None,
            defaultLandmarks=False, defaultTrigPoints=False,
            limitPerLandmark=Parameters.DEFAULT_LIMIT_PER_LANDMARK,
            maxDistance=Parameters.DEFAULT_MAX_DISTANCE,
            minDistance=Parameters.DEFAULT_MIN_DISTANCE,
            distanceBase=Parameters.DEFAULT_DISTANCE_BASE,
            database=None, databaseFile=None,
            databaseFasta=None, subjects=None):
        """
        Use Python function keywords to build an argument parser that can
        used to find or create a database using getDatabaseFromArgs

        @param landmarkNames: a C{list} of C{str} of landmark finder names.
        @param trigPointNames: a C{list} of C{str} of trig finder names.
        @param defaultLandmarks: If C{True}, use the default landmark finders.
        @param defaultTrigPoints: If C{True}, use the default trig point
            finders.
        @param limitPerLandmark: An C{int} limit on the number of pairs to
            yield per landmark.
        @param maxDistance: The C{int} maximum distance permitted between
            yielded pairs.
        @param minDistance: The C{int} minimum distance permitted between
            yielded pairs.
        @param distanceBase: The distance between a landmark and a trig point
            is scaled to be its logarithm using this C{float} base. This
            reduces sensitivity to relatively small differences in distance.
        @param database: An instance of C{light.database.Database}.
        @param databaseFile: The C{str} file name containing a database.
        @param databaseFasta: The name of a FASTA file containing subject
            sequences that should be added to the database.
        @param subjects: A C{dark.reads.Reads} instance containing amino
            acid subject sequences to add to the database.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        parser = argparse.ArgumentParser()
        self.addArgsToParser(parser)
        commandLine = []

        # An in-memory database gets returned immediately, after adding any
        # optional sequences to it.
        if database is not None:
            assert self._allowInMemory, (
                'In-memory database specification not enabled')
            if self._allowPopulation:
                if databaseFasta is not None:
                    for read in FastaReads(databaseFasta, readClass=AARead):
                        database.addSubject(read)
                if subjects is not None:
                    for subject in subjects:
                        database.addSubject(subject)
            return database

        if databaseFile is not None:
            commandLine.extend(['--databaseFile', databaseFile])

        if landmarkNames is not None:
            for landmarkName in landmarkNames:
                commandLine.extend(['--landmark', landmarkName])

        if trigPointNames is not None:
            for trigPointName in trigPointNames:
                commandLine.extend(['--trig', trigPointName])

        if defaultLandmarks:
            commandLine.append('--defaultLandmarks')

        if defaultTrigPoints:
            commandLine.append('--defaultTrigPoints')

        commandLine.extend(['--limitPerLandmark', str(limitPerLandmark),
                            '--maxDistance', str(maxDistance),
                            '--minDistance', str(minDistance),
                            '--distanceBase', str(distanceBase)])

        if self._allowPopulation and databaseFasta is not None:
            commandLine.extend(['--databaseFasta', databaseFasta])

        database = self.getDatabaseFromArgs(parser.parse_args(commandLine))

        if self._allowPopulation and subjects is not None:
            for subject in subjects:
                database.addSubject(subject)

        return database
