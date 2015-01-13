import sys
from json import dumps
from Bio.ExPASy import Prosite


def patternToRegex(pattern):
    """
    Translates the pattern of a prosite entry into a regex.

    @param pattern: the C{str} pattern of a prosite entry.
    """
    switch = {'{': '[^',
              '}': ']',
              'x': '[A-Z]',
              '<': '',
              '>': '',
              '-': '',
              }
    aa = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M',
          'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'X', '[', ']', '(', ')',
          ',', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
    regex = []
    for element in pattern[:-1]:
        if element not in aa:
            regex.append(switch[element])
        else:
            regex.append(element)
    return ''.join(regex)


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
    """
    print 'HERE'
    for record in Prosite.parse(open(prositeDb)):
        print 'InREf', record
        accession = record.accession
        pattern = patternToRegex(record.pattern)
        print 'PATTERN', pattern, accession
        if pattern:
            print >>fp, dumps(
                {
                    'accession': accession,
                    'pattern': pattern,
                }, separators=(',', ':'))

    return fp
