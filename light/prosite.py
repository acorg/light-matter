import sys
from json import dumps
from Bio.ExPASy import Prosite


def patternToRegex(pattern):
    """
    Translates the pattern of a prosite entry into a regex.

    @param pattern: the C{str} pattern of a prosite entry.
    @return: A C{str} that can be used as a regex to match C{pattern}.
    """
    switch = {'{': '[^',
              '}': ']',
              'x': '.',
              '<': '',
              '>': '',
              '-': '',
              }
    aa = set('ARNDCEQGHILKMFPSTWYVUX[](),1234567890')
    return ''.join((element if element in aa else switch[element])
                   for element in pattern)


def prositeToJSON(prositeDb, fp=sys.stdout):
    """
    This is a parser for the prosite database, to turn the relevant entries of
    the prosite database into JSON. Currently it only makes two entries,
    the accession and the pattern. The pattern is converted into a regex.

    The prosite database is available at:
    ftp://ftp.expasy.org/databases/prosite/prosite.dat
    An explanation about the fields and structure of the database is available
    at: http://prosite.expasy.org/prosuser.html

    @param prositeDb: The C{str} filename of the prosite database.
    @param fp: A file pointer.
    @raises AssertionError: if any accession string in the database does not
        start with "PS" or if the database contains a duplicate accession
        string.
    """
    seen = {}
    for record in Prosite.parse(open(prositeDb)):
        accession = record.accession
        assert accession not in seen
        assert accession.startswith('PS')
        seen[accession] = None
        pattern = patternToRegex(record.pattern[:-1])
        if pattern:
            print(dumps(
                {
                    'accession': accession[2:],
                    'pattern': pattern,
                }, separators=(',', ':')), file=fp)
