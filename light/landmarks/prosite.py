import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark
from light.finder import Finder


def loadDatabase():
    filename = join(dirname(light.__file__), '..', 'data',
                    'prosite-20.119.json')
    database = []
    append = database.append
    with open(filename) as fp:
        for line in fp:
            motif = loads(line)
            regex = re.compile(motif['pattern'])
            append({
                'accession': motif['accession'],
                'regex': regex,
            })
    return database

_DATABASE = loadDatabase()


class Prosite(Finder):
    """
    A class for computing statistics based on prosite motifs. The prosite
    database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param dbParams: A C{DatabaseParameters} instance.
    """
    NAME = 'Prosite'
    SYMBOL = 'PS'

    def find(self, read):
        """
        A function that checks if and where a prosite motif in a sequence
        occurs and returns C{Landmark} instances.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for motif in _DATABASE:
            for match in motif['regex'].finditer(read.sequence):
                start = match.start()
                end = match.end()
                length = end - start
                yield Landmark(self.NAME, self.SYMBOL, start, length,
                               motif['accession'])
