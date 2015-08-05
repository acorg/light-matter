import argparse
import sys
from warnings import warn
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from dark.fasta import combineReads, FastaReads
from dark.reads import AARead

from light.connector import SimpleConnector
from light.parameters import Parameters
from light.landmarks import (
    findLandmarks, ALL_LANDMARK_CLASSES, DEFAULT_LANDMARK_CLASSES)
from light.trig import (
    findTrigPoints, ALL_TRIG_CLASSES, DEFAULT_TRIG_CLASSES)
from light.wamp import DEFAULT_URL


class Database:
    """
    An interface to a database connector that maintains a collection of
    sequences (aka subjects) and provide for database insertion and look-up
    operations on them.

    @param params: A C{Parameters} instance.
    @param connector: A C{ligh.SimpleConnector} instance (or other connector)
        for connecting to the database backend(s).
    """
    def __init__(self, params, connector=None):
        self.params = params
        self._connector = connector or SimpleConnector(params)
        # Most of our implementation comes directly from our connector.
        for method in ('addSubject', 'find', 'getIndexBySubject',
                       'getSubjectByIndex', 'getSubjects', 'hashCount',
                       'subjectCount', 'totalResidues', 'totalCoveredResidues',
                       'checksum'):
            setattr(self, method, getattr(self._connector, method))

    def save(self, fp=sys.stdout):
        """
        Save the database state to a file.

        @param fp: A file pointer.
        """
        self.params.save(fp)

        state = {
            '_connectorClassName': self._connector.__class__.__name__,
        }
        dump(state, fp)
        fp.write('\n')
        self._connector.save(fp)

    @classmethod
    def restore(cls, fp=sys.stdin):
        """
        Load a database from a file.

        @param fp: A file pointer.
        @return: An instance of L{Database}.
        @raises ValueError: If a now non-existent connector class name is
            found in the saved database file.
        """
        params = Parameters.restore(fp)
        state = loads(fp.readline()[:-1])

        if state['_connectorClassName'] == SimpleConnector.__name__:
            connector = SimpleConnector.restore(fp)
        else:
            raise ValueError('Unknown backend connector class %r.' %
                             state['_connectorClassName'])

        new = cls(params, connector)

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

        return self._connector.print_(fp, printHashes)


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
        if not (allowCreation or allowFromFile or allowInMemory):
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
        if self._allowFromFile:
            parser.add_argument(
                '--databaseFile',
                help='The name of a file containing a saved database')

        if self._allowWamp:
            parser.add_argument(
                '--debugWamp', action='store_true',
                help='Enable WAMP debug output.')

            parser.add_argument(
                '--wampUrl', default=None,
                help=('The WAMP router URL to connect to. E.g., %s' %
                      DEFAULT_URL))

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
        database = None

        if self._allowFromFile and args.databaseFile:
            with open(args.databaseFile) as fp:
                database = Database.restore(fp)

        elif self._allowWamp and args.wampUrl:
            database = 3  # TODO FIXME

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
