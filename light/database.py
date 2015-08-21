import argparse
import os
import sys
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from Bio.File import as_handle

from dark.fasta import combineReads, FastaReads
from dark.reads import AARead

from light.connector import SimpleConnector, WampServerConnector
from light.backend import Backend
from light.parameters import Parameters
from light.landmarks import (
    findLandmarks, ALL_LANDMARK_CLASSES, DEFAULT_LANDMARK_CLASSES)
from light.trig import (
    findTrigPoints, ALL_TRIG_CLASSES, DEFAULT_TRIG_CLASSES)
from light.wamp import addArgsToParser as addWampArgsToParser


class Database:
    """
    An interface to a database connector that maintains a collection of
    sequences (aka subjects) and provide for database insertion and look-up
    operations on them.

    @param params: A C{Parameters} instance.
    @param connector: A C{ligh.SimpleConnector} instance (or other connector)
        for connecting to the database backend(s).
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """

    SAVE_SUFFIX = '.lmdb'

    def __init__(self, params, connector=None, filePrefix=None):
        self.params = params
        self._connector = connector or SimpleConnector(params)
        self._filePrefix = filePrefix

        # Most of our implementation comes directly from our connector.
        for method in ('addSubject', 'find', 'getIndexBySubject',
                       'getSubjectByIndex', 'getSubjects', 'hashCount',
                       'subjectCount', 'totalResidues', 'totalCoveredResidues',
                       'checksum'):
            setattr(self, method, getattr(self._connector, method))

    def shutdown(self, noSave, filePrefix):
        """
        Shut down the database.

        @param noSave: If C{True}, do not save the database state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        self._connector.shutdown(noSave, filePrefix)

        if not noSave:
            self.save(filePrefix)

    def save(self, fpOrFilePrefix=None):
        """
        Save the database state to a file.

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
            '_connectorClassName': self._connector.__class__.__name__,
        }

        with as_handle(saveFile, 'w') as fp:
            self.params.save(fp)
            dump(state, fp)
            fp.write('\n')

        self._connector.save(fpOrFilePrefix)

    @classmethod
    def restore(cls, fpOrFilePrefix):
        """
        Load a database from a file.

        @param fpOrFilePrefix: A file pointer, or the C{str} prefix of a file
            name, or C{None}. If a C{str}, self.SAVE_SUFFIX is appended to get
            the full file name.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent connector class name is
            found in the saved database file.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + cls.SAVE_SUFFIX
            filePrefix = fpOrFilePrefix
        else:
            saveFile = fpOrFilePrefix
            filePrefix = None

        with as_handle(saveFile) as fp:
            params = Parameters.restore(fp)
            state = loads(fp.readline()[:-1])

        if state['_connectorClassName'] == SimpleConnector.__name__:
            connector = SimpleConnector.restore(fpOrFilePrefix)
        else:
            raise ValueError('Unknown backend connector class %r.' %
                             state['_connectorClassName'])

        new = cls(params, connector, filePrefix=filePrefix)

        return new

    def print_(self, fp=sys.stdout, printHashes=False):
        """
        Print information about the database.

        @param fp: A file pointer to write to.
        @param printHashes: If C{True}, print all hashes and associated
            subjects.
        @return: The value returned by the connector's print_ method.
        """
        self.params.print_(fp)
        print('Connector class: %s' % self._connector.__class__.__name__,
              file=fp)
        print('Subject count: %d' % self.subjectCount(), file=fp)
        print('Hash count: %d' % self.hashCount(), file=fp)
        totalResidues = self.totalResidues()
        print('Total residues: %d' % totalResidues, file=fp)
        if totalResidues:
            print('Coverage: %.2f%%' % (
                100.0 * self.totalCoveredResidues() / totalResidues), file=fp)
        else:
            print('Coverage: 0.00%', file=fp)
        print('Checksum: %d' % self.checksum(), file=fp)

        self._connector.print_(fp, printHashes)


class DatabaseSpecifier:
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
    @param allowWamp: If True, add options that allow the user to specify a
        WAMP-based distributed database.
    @raise ValueError: If the allow options do not permit creation or loading.
    """
    def __init__(self, allowCreation=True, allowPopulation=True,
                 allowInMemory=True, allowFromFile=True, allowWamp=True):
        if not (allowCreation or allowFromFile or allowInMemory or allowWamp):
            raise ValueError('You must either allow database creation, '
                             'loading a database from a file, or passing an '
                             'in-memory database.')
        self._allowCreation = allowCreation
        self._allowPopulation = allowPopulation
        self._allowInMemory = allowInMemory
        self._allowFromFile = allowFromFile
        self._allowWamp = allowWamp

    def addArgsToParser(self, parser):
        """
        Add standard database creation or loading arguments to an argparse
        parser, depending on the allowable ways of specifying a database.

        @param parser: An C{argparse.ArgumentParser} instance.
        """

        parser.add_argument(
            '--filePrefix',
            help=('The prefix of the name of a file containing saved '
                  'data. The suffix "%s" will be added to database '
                  'save files, "%s" to connector save files, and '
                  '"%s-N" to backend save files (where N is a numeric '
                  'backend index).' % (Database.SAVE_SUFFIX,
                                       SimpleConnector.SAVE_SUFFIX,
                                       Backend.SAVE_SUFFIX)))

        if self._allowWamp:
            parser.add_argument(
                '--wamp', action='store_true', default=False,
                help='If True, use a WAMP connected distributed database.')

            addWampArgsToParser(parser)

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
        specify a database: pre-existing in a file (via --filePrefix),
        and then via the creation of a new database (many args). There are
        currently no checks to make sure the user isn't trying to do more
        than one of these at the same time (e.g., using both --filePrefix
        and --defaultLandmarks), the one with the highest priority is silently
        acted on first.

        @param args: Command line arguments as returned by the C{argparse}
            C{parse_args} method.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        def getParametersFromArgs(args):
            landmarkClasses = (
                DEFAULT_LANDMARK_CLASSES if args.defaultLandmarks
                else findLandmarks(args.landmarkFinderNames))

            trigClasses = (
                DEFAULT_TRIG_CLASSES if args.defaultTrigPoints
                else findTrigPoints(args.trigFinderNames))

            if len(landmarkClasses) + len(trigClasses) == 0:
                raise RuntimeError(
                    'Not creating database as no landmark or trig point '
                    'finders were specified. Use --landmark and/or --trig.')

            return Parameters(landmarkClasses, trigClasses,
                              limitPerLandmark=args.limitPerLandmark,
                              maxDistance=args.maxDistance,
                              minDistance=args.minDistance,
                              distanceBase=args.distanceBase)

        database = None

        if self._allowFromFile and args.filePrefix:
            saveFile = args.filePrefix + Database.SAVE_SUFFIX
            if os.path.exists(saveFile):
                database = Database.restore(saveFile)

        elif self._allowWamp and args.wamp:
            params = getParametersFromArgs(args)
            # TODO: what args do we have to pass to the WampServerConnector?
            connector = WampServerConnector(params)
            database = Database(params, connector=connector)

        elif self._allowCreation:
            params = getParametersFromArgs()
            database = Database(params)

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
            database=None, databaseFasta=None, subjects=None, filePrefix=None):
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
        @param databaseFasta: The name of a FASTA file containing subject
            sequences that should be added to the database.
        @param subjects: A C{dark.reads.Reads} instance containing amino
            acid subject sequences to add to the database.
        @param filePrefix: The C{str} prefix of the name of a file containing
            saved data. A suffix will be added to get the various file names
            of the database hash index, the connector, the parameters, etc.
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

        if filePrefix is not None:
            commandLine.extend(['--filePrefix', filePrefix])

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
