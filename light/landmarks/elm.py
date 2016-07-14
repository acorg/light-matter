import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark
from light.finder import Finder


def _loadDatabase():
    """
    Read and convert our JSON version of the eukaryotic linear motif (ELM)
    database.

    The ELM database is available at http://elm.eu.org/downloads

    See data/README.md for details of how to download and convert from the
    original ELM database to our JSON format.

    @return: A C{list} of C{dict}s, each containing 1) an 'identifier' key
        giving the C{str} ELM identifier and 2) a 'pattern' key with a compiled
        regex value.
    """
    filename = join(dirname(light.__file__), '..', 'data',
                    'elm-160713.json')
    database = []
    append = database.append
    with open(filename) as fp:
        for line in fp:
            motif = loads(line)
            pattern = re.compile(motif['pattern'])
            append({
                'identifier': motif['identifier'],
                'pattern': pattern,
            })
    return database


# Read the database (just once) into a singleton for use in the ELM.find
# method.
_DATABASE = _loadDatabase()


class EukaryoticLinearMotif(Finder):
    """
    A class for finding eukariotic linear motifs.
    """
    NAME = 'EukaryoticLinearMotif'
    SYMBOL = 'ELM'

    def find(self, read):
        """
        A function that finds and yields (as C{Landmark}s instances) eukariotic
        linear motifs from a sequence.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for motif in _DATABASE:
            for match in motif['pattern'].finditer(read.sequence):
                start = match.start()
                length = match.end() - start
                yield Landmark(self.NAME, self.SYMBOL, start, length,
                               motif['identifier'])
