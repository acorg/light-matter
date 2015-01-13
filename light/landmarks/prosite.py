import re
from json import loads

from light.features import Landmark


class Prosite(object):
    """
    A class for computing statistics based on prosite motifs. The prosite
    database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param databaseFile: the C{str} name of the prosite database file.
    """
    NAME = 'Prosite'
    SYMBOL = 'P'

    def __init__(self, databaseFile=None):
        if databaseFile:
            self.databaseFile = databaseFile
        else:
            self.databaseFile = '../data/prosite-20-110.json'
        self.database = []
        with open(self.databaseFile) as fp:
            for line in fp:
                self.database.append(loads(line))

    def find(self, read):
        """
        A function that checks if and where a prosite motif in a sequence
        occurs and returns C{Landmark} instances.

        @param read: An instance of C{dark.reads.AARead}.
        @return: A generator that yields C{Landmark} instances.
        """
        for motif in self.database:
            motifRegex = re.compile(motif['pattern'])
            for match in motifRegex.finditer(read.sequence):
                start = match.start()
                end = match.end()
                length = end - start
                yield Landmark(self.NAME, motif['accession'], start,
                               length)
