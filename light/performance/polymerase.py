"""
The Z-scores below are for polymerase sequences taken from Cerny et al., 2014,
'Evolution of Tertiary Structure of Viral RNA Dependent Polymerases.'

Note that the name of 1RW3 has been changed to 4MH8. See
http://www.rcsb.org/pdb/explore/explore.do?structureId=4MH8
In the tables below we use 1RW3
"""

from os.path import join, dirname
from json import loads
from collections import OrderedDict

import light


# The following has long lines and extra embedded spaces to make it simpler
# to read and check for errors. It shoud be identical to Table 2 in Cerny
# et al. (apart from the first line, which is absent in the paper).
_CERNY_TABLE_2 = OrderedDict([
    #           2J7WA 4K6MA 1S49A 1NB6A 3OLBA 1XR7A 1XR6A 3CDWA 2E9ZA 3BSOA 3UQSA 1KHVB 2CKWA 1HI0P 3AXVA 2PUSA 2YI9A 2R7WA 1N35A 3V81C 1MU2A
    ('2J7W-A', []),
    ('4K6M-A', [42.9]),
    ('1S49-A', [22.8, 21.7]),
    ('1NB6-A', [20.5, 17.4, 27.4]),
    ('3OLB-A', [18.1, 16.8, 25.3, 21.5]),
    ('1XR7-A', [18.2, 16.6, 25.1, 20.9, 52.4]),
    ('1XR6-A', [18.0, 16.5, 24.8, 20.7, 52.2, 56.7]),
    ('3CDW-A', [18.0, 16.3, 25.2, 21.0, 53.1, 52.4, 53.1]),
    ('2E9Z-A', [19.2, 17.2, 26.5, 21.6, 41.5, 41.3, 41.0, 41.6]),
    ('3BSO-A', [20.5, 16.5, 27.1, 23.8, 32.0, 32.3, 38.1, 31.8, 32.4]),
    ('3UQS-A', [20.9, 17.7, 28.0, 25.2, 31.1, 31.5, 31.2, 31.4, 32.2, 51.0]),
    ('1KHV-B', [18.7, 17.9, 27.4, 24.3, 32.4, 33.0, 32.9, 33.0, 32.4, 39.3, 42.7]),
    ('2CKW-A', [17.5, 15.0, 24.7, 20.6, 30.4, 30.8, 30.8, 30.9, 30.8, 39.1, 39.4, 43.9]),
    ('1HI0-P', [14.8, 10.6,  4.1, 16.4, 17.2, 17.0, 16.9, 17.7, 15.7, 18.5, 19.1, 17.7, 14.1]),
    ('3AVX-A', [11.1,  7.7, 14.8, 14.1, 14.0, 13.5, 13.6, 14.5, 13.8, 13.2, 14.4, 14.9, 12.6, 12.3]),
    ('2PUS-A', [ 8.4,  6.6, 10.7,  9.5, 12.1, 12.1, 11.9, 12.6, 12.9, 13.4, 13.3, 12.6, 12.9,  9.5,  6.0]),
    ('2YI9-A', [ 9.8,  6.7, 13.9, 12.9, 12.4, 12.3, 12.1, 13.0, 13.5, 15.5, 14.2, 14.0, 13.2, 10.7,  7.7, 42.5]),
    ('2R7W-A', [ 8.9,  9.0, 10.2, 10.5,  9.7,  9.4,  8.3,  8.4,  9.3,  9.4,  9.1, 10.4,  8.5,  9.9,  7.8,  4.6,  4.6]),
    ('1N35-A', [ 6.5,  4.0, 10.3,  7.6,  7.8,  7.3,  7.1,  7.8,  8.1,  7.9,  7.9,  8.1,  8.0,  8.4,  8.0,  6.5,  6.6, 15.4]),
    ('3V81-C', [ 4.7,  1.6,  6.3,  6.5,  5.4,  5.5,  4.9,  4.8,  5.3,  5.5,  5.7,  5.7,  4.9,  3.8,  5.8,  2.8,  2.3,  4.0,  5.9]),
    ('1MU2-A', [ 5.4,  4.0,  7.9,  7.4,  6.2,  6.6,  6.8,  6.9,  6.1,  7.6,  7.9,  6.5,  7.4,  5.5,  7.7,  3.6,  4.3,  4.6,  5.1, 28.5]),
    ('1RW3-A', [ 4.7,  3.4,  7.9,  6.2,  7.2,  7.4,  7.0,  6.8,  6.0,  7.6,  6.8,  7.5,  7.4,  4.9,  6.2,  2.6,  3.0,  4.0,  3.9, 18.2, 20.7]),
    #           2J7WA 4K6MA 1S49A 1NB6A 3OLBA 1XR7A 1XR6A 3CDWA 2E9ZA 3BSOA 3UQSA 1KHVB 2CKWA 1HI0P 3AXVA 2PUSA 2YI9A 2R7WA 1N35A 3V81C 1MU2A
])

# Z_SCORES will be a dict of dicts, indexed first by query id and then by
# subject id, holding the corresponding Z scores (taken from
# _CERNY_TABLE_2).  We'll walk through _CERNY_TABLE_2 to fill in Z_SCORES.
Z_SCORES = {}

for index, queryId in enumerate(_CERNY_TABLE_2):
    Z_SCORES[queryId] = {}
    for subjectId in _CERNY_TABLE_2:
        if queryId != subjectId:
            try:
                score = _CERNY_TABLE_2[subjectId][index]
            except IndexError:
                # The logic here may appear backwards, but it's not. To
                # convince yourself either think about it or look at the
                # testLogic test in test/performance/test_polymerase.py
                score = Z_SCORES[subjectId][queryId]

            Z_SCORES[queryId][subjectId] = score

# Read the BLAST bit score JSON. These values in that file were generated
# by blastp. See the top-level Makefile,
# performance/bin/create-polymerase-json.sh and
# performance/bin/convert-blast-10-to-json.py for details.
_BITSCORES_JSON = join(dirname(dirname(light.__file__)),
                       'performance', 'bit-scores', 'polymerase.json')
BIT_SCORES = loads(open(_BITSCORES_JSON).read())
