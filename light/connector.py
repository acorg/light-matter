import sys

from light.backend import Backend
from light.parameters import Parameters
from light.result import Result


class SimpleConnector:
    """
    Provide a simple connection between a Database and a single Backend.

    @param params: A C{Parameters} instance.
    @param backend: A C{Backend} instance, or C{None} if a backend should be
        created.
    """
    def __init__(self, params, backend=None):
        self.params = params
        self._backend = backend or Backend(params)
        # Most of our implementation comes directly from our backend.
        for method in ('getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects', 'subjectCount', 'hashCount',
                       'totalResidues', 'totalCoveredResidues', 'checksum'):
            setattr(self, method, getattr(self._backend, method))

    def addSubject(self, subject):
        """
        Ask the backend to add a subject.

        @param subject: A C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @return: A tuple of 1) a C{bool} to indicate whether the subject was
            already in the database, and 2) the C{str} subject index.
        """
        preExisting, subjectIndex, _ = self._backend.addSubject(subject)
        return preExisting, subjectIndex

    def find(self, read, significanceMethod=None, scoreMethod=None,
             significanceFraction=None, storeFullAnalysis=False):
        """
        Check which database sequences a read matches.

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

        matches, hashCount, nonMatchingHashes = self._backend.find(
            read, significanceMethod, scoreMethod, significanceFraction,
            storeFullAnalysis)

        return Result(read, self, matches, hashCount,
                      significanceMethod, scoreMethod, significanceFraction,
                      nonMatchingHashes, storeFullAnalysis=storeFullAnalysis)

    def save(self, fp=sys.stdout):
        """
        Save state to a file.

        @param fp: A file pointer.
        @return: The result returned by the 'save' method of our backend.
        """
        self.params.save(fp)
        self._backend.save(fp)

    @classmethod
    def restore(cls, fp=sys.stdin):
        """
        Restore state from a file.

        @param fp: A file pointer.
        @return: An instance of L{SimpleConnector}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved database file. Or if valid JSON cannot be
            loaded from C{fp}.
        """
        return cls(Parameters.restore(fp), Backend.restore(fp))

    def print_(self, fp=sys.stdout, printHashes=False):
        """
        Print information about this connector.

        @param fp: A file pointer to write to.
        @param printHashes: If C{True}, print all hashes and associated
            subjects from the backend.
        """
        if printHashes:
            print('Backends:', file=fp)
            self._backend.print_(fp)
