import six
import argparse
import os
try:
    import asyncio
except ImportError:
    # Python 2
    import trollius as asyncio

try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

import logging

from Bio.File import as_handle

from dark.fasta import combineReads, FastaReads
from dark.reads import AAReadWithX

from light.string import MultilineString
from light.parameters import DatabaseParameters

from light.connector import SimpleConnector
if six.PY3:
    from light.wamp import addArgsToParser as addWampArgsToParser
    from light.connector_wamp import WampServerConnector
    from light.database_wamp import getWampClientDatabase

from light.backend import Backend

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


class Database:
    """
    An interface to a database connector that maintains a collection of
    sequences (aka subjects) and provide for database insertion and look-up
    operations on them.

    @param dbParams: A C{Parameters} instance.
    @param connector: A C{ligh.SimpleConnector} instance (or other connector)
        for connecting to the database backend(s).
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """

    SAVE_SUFFIX = '.lmdb'

    def __init__(self, dbParams, connector=None, filePrefix=None):
        self.dbParams = dbParams
        self._connector = connector or SimpleConnector(dbParams,
                                                       filePrefix=filePrefix)
        self._filePrefix = filePrefix

        # Most of our implementation comes directly from our connector.
        for method in ('addSubject', 'find', 'getIndexBySubject',
                       'getSubjectByIndex', 'getSubjects', 'subjectCount',
                       'totalResidues', 'totalCoveredResidues', 'checksum'):
            setattr(self, method, getattr(self._connector, method))

    def hashCount(self):
        """
        Get the overall hash count for a database.

        @return: An C{int} count of the number of hashes in the database.
        """
        result = self._connector.hashCount()
        if asyncio.iscoroutine(result):
            asyncio.ensure_future(result)
            # TODO: FIX ME!
            raise NotImplementedError()
        else:
            return result

    def shutdown(self, save, filePrefix):
        """
        Shut down the database.

        @param save: If C{True}, save the database state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        self._connector.shutdown(save, filePrefix)

        if save:
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
            self.dbParams.save(fp)
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
            dbParams = DatabaseParameters.restore(fp)
            state = loads(fp.readline()[:-1])

        connectorClassName = state['_connectorClassName']
        if connectorClassName == SimpleConnector.__name__:
            connector = SimpleConnector.restore(fpOrFilePrefix)
        elif six.PY3 and connectorClassName == WampServerConnector.__name__:
            connector = WampServerConnector.restore(fpOrFilePrefix)
        else:
            raise ValueError('Unknown backend connector class %r.' %
                             connectorClassName)

        new = cls(dbParams, connector, filePrefix=filePrefix)

        return new

    def print_(self, printHashes=False, margin=''):
        """
        Print information about the database.

        @param printHashes: If C{True}, print all hashes and associated
            subjects.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the database.
        """
        result = MultilineString(margin=margin)
        append = result.append

        append(self.dbParams.print_(margin=margin), verbatim=True)

        totalResidues = self.totalResidues()
        result.extend([
            'Connector class: %s' % self._connector.__class__.__name__,
            'Subject count: %d' % self.subjectCount(),
            'Hash count: %d' % self.hashCount(),
            'Total residues: %d' % totalResidues,
        ])

        if totalResidues:
            append('Coverage: %.2f%%' % (
                100.0 * self.totalCoveredResidues() / totalResidues))
        else:
            append('Coverage: 0.00%')
        append('Checksum: %d' % self.checksum())

        append('Connector:')
        connector = self._connector.print_(printHashes=printHashes,
                                           margin=margin)
        if connector:
            append(connector, verbatim=True)

        return str(result)


class DatabaseSpecifier:
    """
    Helper class for 1) either creating a new database, loading one from
    saved files, accessing one that's already in memory, or accessing one
    remotely, and then, 2) optionally populating it.

    @param allowCreation: If C{True}, add options that permit the creation
        of a new database.
    @param allowPopulation: If C{True}, add options that allow the automatic
        addition of sequences to the database.
    @param allowInMemory: If C{True}, add the option that allows the user to
        pass an in-memory database (e.g., as constructed in iPython or
        iPythonNotebook or programmatically).
    @param allowWamp: If C{True}, add options that allow the user to specify a
        WAMP-based distributed database. This option can only be used in
        Python 3.
    @raise ValueError: If the allow* options do not permit creation, loading
        from a file, or an in-memory database. Or if the use of WAMP is
        attempted and we're not using Python 3.
    """
    def __init__(self, allowCreation=True, allowPopulation=True,
                 allowInMemory=True, allowWamp=six.PY3):
        if not (allowCreation or allowInMemory or allowWamp):
            raise ValueError('You must either allow database creation, '
                             'loading a database from a file, or passing an '
                             'in-memory database.')

        if allowWamp and not six.PY3:
            raise ValueError('You can only use allowWamp under Python 3.')
        self._allowCreation = allowCreation
        self._allowPopulation = allowPopulation
        self._allowInMemory = allowInMemory
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
                '--wampClient', action='store_true', default=False,
                help=('If True, use a database that is actually just a client '
                      'of a remote WAMP distributed database.'))

            parser.add_argument(
                '--wampServer', action='store_true', default=False,
                help='If True, serve a WAMP-connected distributed database.')

            addWampArgsToParser(parser)

        if self._allowCreation:
            DatabaseParameters.addArgsToParser(parser)

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

    def getDatabaseFromArgs(self, args, dbParams=None):
        """
        Read an existing database (if args.database is given) or create
        one from args.

        There is an order of preference in examining the arguments used to
        specify a database: pre-existing in a file (via --filePrefix),
        and then via the creation of a new database (many args). There are
        currently no checks to make sure the user isn't trying to do
        conflicting things, such as restoring from a file and also specifying
        landmark finders, the one with the highest priority is silently
        acted on first.

        @param args: Command line arguments as returned by the C{argparse}
            C{parse_args} method.
        @param dbParams: A C{DatabaseParameters} instance or C{None} to use
            default parameters.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """

        dbParams = dbParams or DatabaseParameters()
        database = None
        filePrefix = args.filePrefix

        if filePrefix:
            # Check to see which save files exist so we know if a restore
            # is possible, and what kind of database & connector were
            # involved.
            exists = os.path.exists
            dbSaveFile = filePrefix + Database.SAVE_SUFFIX
            scSaveFile = filePrefix + SimpleConnector.SAVE_SUFFIX
            if six.PY3:
                wcSaveFile = filePrefix + WampServerConnector.SAVE_SUFFIX
            beSaveFile = filePrefix + Backend.SAVE_SUFFIX
            if exists(dbSaveFile):
                if exists(scSaveFile):
                    if exists(beSaveFile):
                        database = Database.restore(filePrefix)
                    else:
                        raise RuntimeError(
                            'A database save file (%s) and simple connector '
                            'save file (%s) are both present, but no backend '
                            'save file (%s) exists!' %
                            (dbSaveFile, scSaveFile, beSaveFile))
                elif six.PY3 and exists(wcSaveFile):
                    # We do not have to check for backend save files in the
                    # case of a WAMP database, as these may be on another
                    # machine.
                    database = Database.restore(filePrefix)
                else:
                    if six.PY3:
                        raise RuntimeError(
                            'A database save file (%s) is present, but no '
                            'simple connector save file (%s) or WAMP '
                            'connector save file (%s) exists!' %
                            (dbSaveFile, scSaveFile, wcSaveFile))
                    else:
                        raise RuntimeError(
                            'A database save file (%s) is present, but no '
                            'simple connector save file (%s) exists!' %
                            (dbSaveFile, scSaveFile))

        if database is None and self._allowWamp:
            if args.wampServer:
                dbParams = DatabaseParameters.fromArgs(args)
                connector = WampServerConnector(dbParams,
                                                filePrefix=filePrefix)
                database = Database(dbParams, connector=connector,
                                    filePrefix=filePrefix)
            elif args.wampClient:
                database = getWampClientDatabase(args)

        if database is None and self._allowCreation:
            # A new in-memory database, with a simple connector and a local
            # backend.
            dbParams = DatabaseParameters.fromArgs(args)
            database = Database(dbParams, filePrefix=filePrefix)

        if database is None and self._allowWamp:
            # Last try: guess that they might be wanting to talk to a WAMP
            # database, even though --wampClient isn't specified.
            database = getWampClientDatabase(args)

        if database is None:
            raise RuntimeError(
                'Not enough information given to specify a database, and no '
                'remote WAMP database could be found.')

        if self._allowPopulation:
            for read in combineReads(
                    args.databaseFasta, args.databaseSequences,
                    readClass=AAReadWithX, upperCase=True):
                database.addSubject(read)

        return database

    def getDatabaseFromKeywords(self, database=None, databaseFasta=None,
                                subjects=None, filePrefix=None,
                                wampServer=False, wampClient=False, **kwargs):
        """
        Use Python function keywords to build an argument parser that can
        used to find or create a database using getDatabaseFromArgs

        @param database: An instance of C{light.database.Database}.
        @param databaseFasta: The name of a FASTA file containing subject
            sequences that should be added to the database.
        @param subjects: A C{dark.reads.Reads} instance containing amino
            acid subject sequences to add to the database.
        @param filePrefix: The C{str} prefix of the name of a file containing
            saved data. A suffix will be added to get the various file names
            of the database hash index, the connector, the parameters, etc.
        @param wampServer: If C{True} the database should act as WAMP server.
        @param wampClient: If C{True} create a database that is actually a
            client connection to a WAMP database.
        @param kwargs: Additional arguments specifying values for database
            parameters (some of which may be used by specific feature finders).
            To see the available database parameters, see the docstring of
            the DatabaseParameters class.
        @raise ValueError: If a database cannot be found or created.
        @return: A C{light.database.Database} instance.
        """
        # An pre-existing in-memory database gets returned immediately,
        # after adding any optional sequences to it.
        if database is not None:
            assert self._allowInMemory, (
                'In-memory database specification not enabled.')
        else:
            # Use self.getDatabaseFromArgs to get the database.
            parser = argparse.ArgumentParser()
            self.addArgsToParser(parser)

            dbParams = DatabaseParameters(**kwargs)
            print('landmarkFinders:', dbParams.landmarkFinders)
            print('trigPointFinders:', dbParams.trigPointFinders)
            commandLine = []

            if filePrefix is not None:
                commandLine.extend(['--filePrefix', filePrefix])

            if wampServer:
                commandLine.append('--wampServer')

            if wampClient:
                commandLine.append('--wampClient')

            args = parser.parse_args(commandLine)
            database = self.getDatabaseFromArgs(args, dbParams)

        if self._allowPopulation:
            if databaseFasta:
                for read in FastaReads(databaseFasta,
                                       readClass=AAReadWithX,
                                       upperCase=True):
                    database.addSubject(read)
            if subjects is not None:
                for subject in subjects:
                    database.addSubject(subject)

        return database
