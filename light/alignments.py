import bz2
from json import loads
import copy

from dark.reads import AARead, Reads
from dark.hsp import HSP
from dark.alignments import (
    Alignment, ReadsAlignments, ReadsAlignmentsParams, ReadAlignments)
from dark.utils import numericallySortFilenames

from light.params import checkCompatibleParams


class LightAlignment(Alignment):
    """
    Hold information about a light matter read alignment. This is the same as
    a dark matter alignment (of which it is a subclass), with the addition of
    match details for landmarks and trig points.

    @param subjectLength: The C{int} length of the sequence a read matched
        against.
    @param subjectTitle: The C{str} title of the sequence a read matched
        against.
    """


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


def jsonDictToAlignments(lightDict, database):
    """
    Take a dict of light matter results (converted from JSON) for a single
    read and produce a list of alignments.

    @param lightDict: A C{dict}, created from a line of light matter result
        JSON.
    @param database: A C{light.database.Database} instance.
    @return: A C{list} of L{light.alignment.LightAlignment} instances.
    """
    lightAlignments = []

    for alignment in lightDict['alignments']:
        subjectId, subjectSequence = database.subjectInfo[
            alignment['subjectIndex']]
        lightAlignment = LightAlignment(len(subjectSequence), subjectId)
        for hsp in alignment['hsps']:
            h = HSP(hsp['score'])
            # Ugh, manually put the additional light matter HSP info onto
            # the HSP instance. For now.
            h.hspInfo = hsp['hspInfo']
            lightAlignment.addHsp(h)
        lightAlignments.append(lightAlignment)

    return lightAlignments


class JSONRecordsReader(object):
    """
    Provide a method that yields JSON records from a file. Store and
    make accessible the run-time parameters.

    NOTE: this class currently has no tests!

    @param filename: A C{str} filename containing JSON light matter records.
    @param database: A C{light.database.Database} instance.
    """
    def __init__(self, filename, database):
        self._filename = filename
        self._database = database
        self._open(filename)

    def _open(self, filename):
        """
        Open an input file. Set self._fp to point to it. Read the first
        line of parameters.

        @param filename: A C{str} filename containing JSON light matter
            results.
        @raise ValueError: if the first line of the file isn't valid JSON,
            if the input file is empty or if the database and output file
            checksums don't match.
        """
        if filename.endswith('.bz2'):
            self._fp = bz2.BZ2File(filename)
        else:
            self._fp = open(filename)

        line = self._fp.readline()
        if not line:
            raise ValueError('JSON file %r was empty.' % self._filename)

        try:
            self.params = loads(line[:-1])
        except ValueError as e:
            raise ValueError(
                'Could not convert first line of %r to JSON (%s). '
                'Line is %r.' % (self._filename, e, line[:-1]))

        if self.params['checksum'] != self._database.checksum:
            raise ValueError(
                'Database and output file have different checksums.')

    def readAlignments(self):
        """
        Read lines of JSON from self._filename, convert them to read alignments
        and yield them.

        @raise ValueError: If any of the lines in the file cannot be converted
            to JSON.
        @return: A generator that yields C{dark.alignments.ReadAlignments}
            instances.
        """
        if self._fp is None:
            self._open(self._filename)

        try:
            for lineNumber, line in enumerate(self._fp, start=2):
                try:
                    record = loads(line[:-1])
                except ValueError as e:
                    raise ValueError(
                        'Could not convert line %d of %r to JSON (%s). '
                        'Line is %r.' %
                        (lineNumber, self._filename, e, line[:-1]))
                else:
                    read = AARead(record['queryId'], record['querySequence'])
                    alignments = jsonDictToAlignments(record, self._database)
                    yield ReadAlignments(read, alignments)
        finally:
            self._fp.close()
            self._fp = None
