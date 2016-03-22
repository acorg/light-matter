import logging

from Bio.File import as_handle

from light.backend import Backend
from light.parameters import DatabaseParameters
from light.result import Result
from light.string import MultilineString

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


class SimpleConnector(object):
    """
    Provide a simple in-memory connection between a Database and a single
    Backend.

    @param dbParams: A C{DatabaseParameters} instance.
    @param backend: A C{Backend} instance, or C{None} if a backend should be
        created.
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """

    SAVE_SUFFIX = '.lmco'

    def __init__(self, dbParams, backend=None, filePrefix=None):
        self.dbParams = dbParams
        if backend:
            self._backend = backend
        else:
            self._backend = Backend(filePrefix=filePrefix)
            self._backend.configure(dbParams)
        self._filePrefix = filePrefix

        # Most of our implementation comes directly from our backend.
        for method in ('addSubject', 'getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects', 'subjectCount', 'hashCount',
                       'totalResidues', 'totalCoveredResidues', 'checksum'):
            setattr(self, method, getattr(self._backend, method))

    def find(self, read, findParams=None, storeFullAnalysis=False,
             subjectIndices=None):
        """
        Check which database sequences a read matches.

        @param read: A C{dark.read.AARead} instance.
        @param findParams: An instance of C{light.parameters.FindParameters} or
            C{None} to use default find parameters.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @param subjectIndices: A C{set} of subject indices, or C{None}. If a
            set is passed, only subject indices in the set will be returned
            in the results. If C{None}, all matching subject indices are
            returned.
        @return: A C{light.result.Result} instance.
        """
        matches, hashCount, nonMatchingHashes = self._backend.find(
            read, storeFullAnalysis=storeFullAnalysis,
            subjectIndices=subjectIndices)

        return Result(read, self, matches, hashCount, findParams,
                      nonMatchingHashes=nonMatchingHashes,
                      storeFullAnalysis=storeFullAnalysis)

    def shutdown(self, save, filePrefix):
        """
        Shut down the connector.

        @param save: If C{True}, save the connector state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        self._backend.shutdown(save, filePrefix)

        if save:
            self.save(filePrefix)

    def save(self, fpOrFilePrefix=None):
        """
        Save state to a file.

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

        with as_handle(saveFile, 'w') as fp:
            self.dbParams.save(fp)

        self._backend.save(fpOrFilePrefix)

    @classmethod
    def restore(cls, fpOrFilePrefix):
        """
        Restore state from a file.

        @param fpOrFilePrefix: A file pointer or the C{str} prefix of a file
            name. If a C{str}, self.SAVE_SUFFIX is appended to get the full
            file name.
        @return: An instance of L{SimpleConnector}.
        @raises ValueError: If valid JSON cannot be loaded from C{fp}.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + cls.SAVE_SUFFIX
            filePrefix = fpOrFilePrefix
        else:
            saveFile = fpOrFilePrefix
            filePrefix = None

        with as_handle(saveFile) as fp:
            dbParams = DatabaseParameters.restore(fp)

        return cls(dbParams, backend=Backend.restore(fpOrFilePrefix),
                   filePrefix=filePrefix)

    def print_(self, printHashes=False, margin='', result=None):
        """
        Print information about this connector.

        @param printHashes: If C{True}, print all hashes and associated
            subjects from the backend.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} representation of
            the connector, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        if printHashes:
            result.append('Backends:')
            result.indent()
            self._backend.print_(margin=margin, result=result)
            result.outdent()

        if not returnNone:
            return str(result)
