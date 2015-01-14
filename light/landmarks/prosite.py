import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark


class Prosite(object):
    """
    A class for computing statistics based on prosite motifs. The prosite
    database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param databaseFile: the C{str} name of the prosite database file.
    @raises AssertionError: if any accession string in the database does not
        start with "PS".
    """
    NAME = 'Prosite'
    SYMBOL = 'PS'

    def __init__(self, databaseFile=None):
        dbFile = databaseFile or join(dirname(light.__file__), '..', 'data',
                                      'prosite-20-110.json')
        self.database = []
        seen = {}
        with open(dbFile) as fp:
            for line in fp:
                motif = loads(line)
                regex = re.compile(motif['pattern'])
                assert motif['accession'] not in seen
                seen[motif['accession']] = None
                assert motif['accession'].startswith('PS')
                self.database.append({
                                     'accession': motif['accession'][2:],
                                     'regex': regex,
                                     })

    def find(self, read):
        """
        A function that checks if and where a prosite motif in a sequence
        occurs and returns C{Landmark} instances.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        @raise: If the first two letters of the accession are not 'PS', raise.
        """
        for motif in self.database:
            for match in motif['regex'].finditer(read.sequence):
                start = match.start()
                end = match.end()
                length = end - start
                yield Landmark(self.NAME, self.SYMBOL, start,
                               length, motif['accession'])
