import bz2
from json import loads

from dark.reads import AARead
from dark.hsp import HSP
from dark.alignments import ReadAlignments

from light.alignments import LightAlignment
# from light.hsp import normalizeHSP


class JSONRecordsReader(object):
    """
    Provide a method that yields JSON records from a file. Store, check, and
    make accessible the run-time parameters.

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
        @raise ValueError: if the first line of the file isn't valid JSON
            or if the input file is empty.
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

    def _dictToAlignments(self, lightDict):
        """
        Take a dict of light matter results for a single read and convert it
        to a list of alignments.

        @param lightDict: A C{dict}, created from a line of light result JSON.
        @return: A C{list} of L{light.alignment.LightAlignment} instances.
        """
        lightAlignments = []

        for alignment in lightDict['alignments']:
            subjectId, subjectSequence = self._database.subjectInfo[
                alignment['subjectIndex']]
            lightAlignment = LightAlignment(len(subjectSequence), subjectId,
                                            alignment['matchInfo'])
            lightAlignment.addHsp(HSP(alignment['matchScore']))
            lightAlignments.append(lightAlignment)

        return lightAlignments

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
                    alignments = self._dictToAlignments(record)
                    yield ReadAlignments(read, alignments)
        finally:
            self._fp.close()
            self._fp = None
