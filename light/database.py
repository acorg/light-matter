import argparse
import sys
from warnings import warn
from collections import defaultdict
from ujson import dump, dumps, load
from binascii import crc32
from operator import attrgetter

from dark.fasta import combineReads

from light.reads import ScannedRead
from light.result import Result
from light.landmarks import (
    findLandmark, findLandmarks, ALL_LANDMARK_CLASSES,
    DEFAULT_LANDMARK_CLASSES)
from light.trig import (
    findTrigPoint, findTrigPoints, ALL_TRIG_CLASSES, DEFAULT_TRIG_CLASSES)


class Database(object):
    """
    Maintain a collection of sequences ("subjects") and provide for database
    (index, search) operations on them.

    @param landmarkClasses: A C{list} of landmark finder classes, or C{None}
        to use the default set of landmark finder classes.
    @param trigPointClasses: A C{list} of trig point finder classes, or C{None}
        to use the default set of trig point finder classes.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param bucketFactor: A C{int} factor by which the distance between
        landmark and trig point is divided, to influence sensitivity.
    """

    # Database construction and look-up defaults. See explanations in
    # docstring above.
    DEFAULT_LIMIT_PER_LANDMARK = 10
    DEFAULT_MAX_DISTANCE = 200
    DEFAULT_MIN_DISTANCE = 1
    DEFAULT_BUCKET_FACTOR = 1  # Should be >1. Can we use a float?

    # The default fraction of all (landmark, trig point) pairs for a
    # scannedRead that need to fall into the same offset delta histogram
    # bucket for that bucket to be considered a significant match with a
    # database title.
    DEFAULT_SIGNIFICANCE_FRACTION = 0.25

    def __init__(self, landmarkClasses, trigPointClasses,
                 limitPerLandmark=None, maxDistance=None, minDistance=None,
                 bucketFactor=None):
        self.landmarkClasses = (
            DEFAULT_LANDMARK_CLASSES if landmarkClasses is None
            else landmarkClasses)

        self.trigPointClasses = (
            DEFAULT_TRIG_CLASSES if trigPointClasses is None
            else trigPointClasses)

        self.limitPerLandmark = (
            self.DEFAULT_LIMIT_PER_LANDMARK if limitPerLandmark is None
            else limitPerLandmark)

        self.maxDistance = (
            self.DEFAULT_MAX_DISTANCE if maxDistance is None
            else maxDistance)

        self.minDistance = (
            self.DEFAULT_MIN_DISTANCE if minDistance is None
            else minDistance)

        self.bucketFactor = (
            self.DEFAULT_BUCKET_FACTOR if bucketFactor is None
            else bucketFactor)

        if self.bucketFactor <= 0:
            raise ValueError('bucketFactor must be > 0.')

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
        for landmarkClass in self.landmarkClasses:
            self.landmarkFinders.append(landmarkClass())
        self.trigPointFinders = []
        for trigPointClass in self.trigPointClasses:
            self.trigPointFinders.append(trigPointClass())
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
        distance = ((trigPoint.offset - landmark.offset)
                    // self.bucketFactor)
        return '%s:%s:%s' % (landmark.hashkey(), trigPoint.hashkey(),
                             distance)

    def scan(self, sequence):
        """
        Makes an instance of C{light.reads.ScannedRead}.

        @param sequence: a C{dark.read.AARead} instance.
        @return: a C{light.reads.ScannedRead} instance.
        """
        scannedSequence = ScannedRead(sequence)
        for landmarkFinder in self.landmarkFinders:
            for landmark in landmarkFinder.find(sequence):
                scannedSequence.landmarks.append(landmark)

        for trigFinder in self.trigPointFinders:
            for trigPoint in trigFinder.find(sequence):
                scannedSequence.trigPoints.append(trigPoint)
        return scannedSequence

    def getScannedPairs(self, scannedSequence):
        """
        Get the (landmark, trigPoint) pairs from a ScannedRead instance.

        @param scannedSequence: A C{light.reads.ScannedRead} instance.
        @return: C{light.reads.ScannedRead.getPairs}
        """
        return scannedSequence.getPairs(limitPerLandmark=self.limitPerLandmark,
                                        maxDistance=self.maxDistance,
                                        minDistance=self.minDistance)

    def addSubject(self, subject):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: a C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @return: The C{int} subject index of the added subject.
        """
        subjectInfo = (subject.id, subject.sequence)
        self._updateChecksum(subjectInfo)
        self.subjectInfo.append(subjectInfo)
        subjectIndex = self.subjectCount
        self.subjectCount += 1
        self.totalResidues += len(subject)

        scannedSubject = self.scan(subject)

        self.totalCoveredResidues += len(scannedSubject.coveredIndices())

        for landmark, trigPoint in self.getScannedPairs(scannedSubject):
            key = self.key(landmark, trigPoint)
            try:
                subjectDict = self.d[key]
            except KeyError:
                self.d[key] = subjectDict = {}

            try:
                subjectDict[str(subjectIndex)].append(landmark.offset)
            except KeyError:
                subjectDict[str(subjectIndex)] = [landmark.offset]

        return subjectIndex

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
            'landmarkClasses': [cls.NAME for cls in self.landmarkClasses],
            'trigPointClasses': [cls.NAME for cls in self.trigPointClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'minDistance': self.minDistance,
            'subjectCount': self.subjectCount,
            'totalResidues': self.totalResidues,
            'totalCoveredResidues': self.totalCoveredResidues,
            'bucketFactor': self.bucketFactor,
        })

        return fp

    def find(self, read, significanceFraction=None, storeFullAnalysis=False):
        """
        A function which takes a read, computes all hashes for it, looks up
        matching hashes and checks which database sequence it matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if significanceFraction is None:
            significanceFraction = self.DEFAULT_SIGNIFICANCE_FRACTION

        scannedRead = self.scan(read)

        matches = defaultdict(list)
        nonMatchingHashes = set()
        hashCount = 0

        for landmark, trigPoint in self.getScannedPairs(scannedRead):
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
                if storeFullAnalysis:
                    nonMatchingHashes.add((landmark, trigPoint))
            else:
                for subjectIndex, subjectOffsets in subjectDict.iteritems():
                    subjectIndex = int(subjectIndex)
                    subjectLength = len(self.subjectInfo[subjectIndex][1])
                    matches[subjectIndex].append({
                        'landmark': landmark,
                        'subjectLength': subjectLength,
                        'subjectOffsets': subjectOffsets,
                        'trigPoint': trigPoint,
                    })

        return Result(scannedRead, matches, hashCount, significanceFraction,
                      self.bucketFactor, nonMatchingHashes,
                      storeFullAnalysis=storeFullAnalysis)

    def save(self, fp=sys.stdout):
        """
        Save the database to a file.

        @param fp: A file pointer.
        """
        state = {
            'checksum': self.checksum,
            'landmarkClassNames': [cls.NAME for cls in self.landmarkClasses],
            'trigPointClassNames': [cls.NAME for cls in self.trigPointClasses],
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

        landmarkClasses = []
        for landmarkClassName in state['landmarkClassNames']:
            cls = findLandmark(landmarkClassName)
            if cls:
                landmarkClasses.append(cls)
            else:
                print >>sys.stderr, (
                    'Could not find landscape finder class %r! Has that '
                    'class been renamed or removed?' % landmarkClassName)
                sys.exit(1)

        trigPointClasses = []
        for trigPointClassName in state['trigPointClassNames']:
            cls = findTrigPoint(trigPointClassName)
            if cls:
                trigPointClasses.append(cls)
            else:
                print >>sys.stderr, (
                    'Could not find trig point finder class %r! Has that '
                    'class been renamed or removed?' % trigPointClassName)
                sys.exit(1)

        database = Database(landmarkClasses, trigPointClasses,
                            limitPerLandmark=state['limitPerLandmark'],
                            maxDistance=state['maxDistance'],
                            minDistance=state['minDistance'],
                            bucketFactor=state['bucketFactor'])

        # Monkey-patch the new database instance to restore its state.
        for attr in ('checksum', 'd', 'subjectCount', 'totalResidues',
                     'totalCoveredResidues', 'subjectInfo'):
            setattr(database, attr, state[attr])

        return database


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
            raise ValueError('You must either allow database creation, or '
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
        if self._allowInMemory:
            parser.add_argument(
                '--database',
                help=('An in-memory database (e.g., as created in ipython or '
                      'iPythonNotebook'))

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
                default=Database.DEFAULT_LIMIT_PER_LANDMARK,
                help=('A limit on the number of pairs to yield per landmark '
                      'per read.'))

            parser.add_argument(
                '--maxDistance', type=int,
                default=Database.DEFAULT_MAX_DISTANCE,
                help='The maximum distance permitted between yielded pairs.')

            parser.add_argument(
                '--minDistance', type=int,
                default=Database.DEFAULT_MIN_DISTANCE,
                help='The minimum distance permitted between yielded pairs.')

            parser.add_argument(
                '--bucketFactor', type=int,
                default=Database.DEFAULT_BUCKET_FACTOR,
                help=('A factor by which the distance between landmark and '
                      'trig point is divided.'))

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
        specify a database: pre-existing in memory (via --database),
        pre-existing in a file (via --databaseFile), and lastly via the
        creation of a new database (many args). There are currently no checks
        to make sure the user isn't trying to do more than one of these at
        the same time (e.g., using both --database and --databaseFile), the
        one with the highest priority is silently acted on first.

        @param args: Command line arguments as returned by the C{argparse}
            C{parse_args} method.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        if self._allowInMemory and args.database:
            database = args.database

        elif self._allowFromFile and args.databaseFile:
            with open(args.databaseFile) as fp:
                database = Database.load(fp)

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

            database = Database(landmarkClasses, trigClasses,
                                args.limitPerLandmark, args.maxDistance,
                                args.minDistance, args.bucketFactor)

        if self._allowPopulation:
            for read in combineReads(args.databaseFasta,
                                     args.databaseSequences):
                database.addSubject(read)

        return database

    def getDatabaseFromKeywords(
            self, landmarkNames=None, trigPointNames=None,
            defaultLandmarks=False, defaultTrigPoints=False,
            limitPerLandmark=Database.DEFAULT_LIMIT_PER_LANDMARK,
            maxDistance=Database.DEFAULT_MAX_DISTANCE,
            minDistance=Database.DEFAULT_MIN_DISTANCE,
            bucketFactor=Database.DEFAULT_BUCKET_FACTOR,
            database=None, databaseFile=None,
            databaseFasta=None, databaseSequences=None):
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
        @param bucketFactor: A C{int} factor by which the distance between
            landmark and trig point is divided.
        @param database: An instance of C{light.database.Database}.
        @param databaseFile: The C{str} file name containing a database.
        @param databaseFasta: The name of a FASTA file containing the
            sequences that should be added to the database.
        @param databaseSequences: A C{list} of amino acid sequences to add
            to the database. The sequence id will be the text up to the last
            space, if any, otherwise will be automatically assigned.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        parser = argparse.ArgumentParser()
        self.addArgsToParser(parser)
        commandLine = []

        if database is not None:
            commandLine.extend(['--database', database])

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
                            '--bucketFactor', str(bucketFactor)])

        if databaseFasta is not None:
            commandLine.extend(['--databaseFasta', databaseFasta])

        if databaseSequences is not None:
            for databaseSequence in databaseSequences:
                commandLine.extend(['--databaseSequence', databaseSequence])

        return self.getDatabaseFromArgs(parser.parse_args(commandLine))
