import copy

from dark.alignments import ReadsAlignments, ReadsAlignmentsParams
from dark.utils import numericallySortFilenames
from dark.reads import AARead, Reads

from light.params import checkCompatibleParams
from light.conversion import JSONRecordsReader


class LightReadsAlignments(ReadsAlignments):
    """
    Hold information about a set of light matter results.

    @param resultFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing our light matter output
        (possibly bzip2 compressed).
    @param database: A C{light.database.Database} instance.
    @param sortFilenames: A C{bool}. If C{True}, C{resultFilenames} will be
        sorted by numeric prefix (using L{numericallySortFilenames}) before
        being read. This can be used to conveniently sort the files produced
        by our HTCondor jobs.
    @raises ValueError: if a file type is not recognized, or if light matter
        parameters in all files do not match.
    """

    def __init__(self, resultFilenames, database, sortFilenames=True):
        if type(resultFilenames) == str:
            resultFilenames = [resultFilenames]
        if sortFilenames:
            self.resultFilenames = numericallySortFilenames(resultFilenames)
        else:
            self.resultFilenames = resultFilenames
        self._database = database

        # Add a dictionary that will allow us to look up databae subjects
        # by title.
        self._subjects = dict(database.subjectInfo)

        # Prepare application parameters in order to initialize self.
        self._reader = self._getReader(self.resultFilenames[0])
        params = copy.deepcopy(self._reader.params)

        applicationParams = ReadsAlignmentsParams(
            'light', applicationParams=params, subjectIsNucleotides=False)

        # We will add to self._reads as we go through the results (which
        # contain the read ids and sequences). Note that this means the
        # reads will not be available to the ReadsAlignments instance until
        # after all results have been read.
        self._reads = Reads()

        ReadsAlignments.__init__(self, self._reads, applicationParams)

    def _getReader(self, filename):
        """
        Obtain a JSON record reader for light matter results.

        @param filename: The C{str} file name holding the JSON.
        """
        if filename.endswith('.json') or filename.endswith('.json.bz2'):
            return JSONRecordsReader(filename, self._database)
        else:
            raise ValueError(
                'Unknown light matter result file suffix for file %r.' %
                filename)

    def iter(self):
        """
        Extract light matter results and yield C{ReadAlignments} instances.

        For each file except the first, check that the parameters are
        compatible with those found (above, in __init__) in the first file.

        @return: A generator that yields C{ReadAlignments} instances.
        """
        # Note that self._reader is already initialized (in __init__) for
        # the first input file. This is less clean than it could be, but it
        # makes testing easier, since open() is then only called once for
        # each input file.

        reader = self._reader
        first = True

        for resultFilename in self.resultFilenames:
            if first:
                # No need to check params in the first file. We already read
                # them in and stored them in __init__.
                first = False
            else:
                reader = self._getReader(resultFilename)
                differences = checkCompatibleParams(
                    self.params.applicationParams, reader.params)
                if differences:
                    raise ValueError(
                        'Incompatible light matter parameters found. The '
                        'parameters in %s differ from those originally found '
                        'in %s. %s' %
                        (resultFilename, self.resultFilenames[0], differences))

            for readAlignments in reader.readAlignments():
                yield readAlignments

    def getSubjectSequence(self, title):
        """
        Obtain information about a subject sequence, given its title.

        @param title: A C{str} sequence title from a light matter match.
        @return: An C{AARead} instance.
        """
        return AARead(title, self._subjects[title])
