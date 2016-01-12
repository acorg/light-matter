import re
from json import loads
from os.path import dirname, join

import light
from light.features import Landmark, Finder


_DB_FILE = join(dirname(light.__file__), '..', 'data', 'prosite-20.119.json')


class Prosite(Finder):
    """
    A class for computing statistics based on prosite motifs. The prosite
    database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param featureLengthBase: See the light.features.Finder class for
        documentation of this parameter. It is unused in the Prosite finder
        but we need it in our signature because all finders must be able to
        accept a featureLengthBase argument.
    """
    NAME = 'Prosite'
    SYMBOL = 'PS'

    def __init__(self, featureLengthBase=1.35):
        Finder.__init__(self, featureLengthBase)
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
