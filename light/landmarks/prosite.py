import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark, Finder


_DB_FILE = join(dirname(light.__file__), '..', 'data', 'prosite-20-110.json')


class Prosite(Finder):
    """
    A class for computing statistics based on prosite motifs. The prosite
    database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param distanceBase: A C{float} by which the length of the landmark or
        trigPoint will be divided. The default of 1.0 set below is only used by
        the tests. Otherwise, it is set by the light.database.Database class.
    """
    NAME = 'Prosite'
    SYMBOL = 'PS'

    def __init__(self, distanceBase=1.0):
        Finder.__init__(self, distanceBase)
        self.database = []
        with open(_DB_FILE) as fp:
            for line in fp:
                motif = loads(line)
                regex = re.compile(motif['pattern'])
                self.database.append({
                                     'accession': motif['accession'],
                                     'regex': regex,
                                     })

    def find(self, read):
        """
        A function that checks if and where a prosite motif in a sequence
        occurs and returns C{Landmark} instances.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for motif in self.database:
            for match in motif['regex'].finditer(read.sequence):
                start = match.start()
                end = match.end()
                length = end - start
                yield Landmark(self.NAME, self.SYMBOL, start, length,
                               motif['accession'])
