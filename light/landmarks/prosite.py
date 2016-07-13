import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark
from light.finder import Finder


def _loadDatabase():
    """
    Read and convert our JSON version of the Prosite database.

    The prosite database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat An explanation about
    the fields and structure of the database is available at:
    http://prosite.expasy.org/prosuser.html

    See data/README.md for details of automated download and conversion from
    the Prosite file to our JSON format.

    @return: A C{list} of C{dict}s, each containing 1) an 'accession' key
        giving the C{str} Prosite accession number and 2) a 'regex' key with
        a compiled regex value.
    """
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


# Read the database (just once) into a singleton for use in the
# Prosite.find method.
_DATABASE = _loadDatabase()


class Prosite(Finder):
    """
    A class for finding prosite motifs.
    """
    NAME = 'Prosite'
    SYMBOL = 'PS'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) Prosite
        motifs from a sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for motif in _DATABASE:
            for match in motif['regex'].finditer(read.sequence):
                start = match.start()
                length = match.end() - start
                yield Landmark(self.NAME, self.SYMBOL, start, length,
                               motif['accession'])
