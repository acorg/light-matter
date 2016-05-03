"""
The Z-scores below are for ha sequences.
"""

from collections import OrderedDict


# The following has long lines and extra embedded spaces to make it simpler
# to read and check for errors.
_Z_SCORES_TABLE = OrderedDict([
    #                     4f3ze 1ru7a 1ruyh 3ztna 2wr7c 5hmge 2fk0k 4dj6a 1ti8a 1jsda 4i78b 4k3xe 4fqkc
    ('pdb_4f3z_e_H1N2',   []),
    ('pdb_1ru7_a_H1N1',   [37.7]),
    ('pdb_1ruy_h_H1N1',   [37.9, 42.0]),
    ('pdb_3ztn_a_H1N1',   [36.6, 38.1, 38.4]),
    ('pdb_2wr7_c_H2N2',   [33.0, 34.7, 34.1, 35.9]),
    ('pdb_5hmg_e_H3N2',   [29.9, 31.9, 31.5, 30.4, 28.1]),
    ('pdb_2fk0_k_H5N1',   [35.6, 36.6, 36.5, 37.9, 36.6, 29.9]),
    ('pdb_4dj6_a_H7N7',   [30.1, 32.0, 31.4, 30.6, 28.8, 31.7, 30.2]),
    ('pdb_1ti8_a_H7N3',   [29.7, 31.6, 31.0, 30.1, 28.2, 31.4, 29.9, 40.4]),
    ('pdb_1jsd_a_H9N2',   [32.7, 35.8, 35.2, 33.5, 31.0, 30.6, 32.5, 30.9, 30.2]),
    ('pdb_4i78_b_H17N10', [31.8, 33.4, 32.9, 34.3, 32.3, 27.6, 33.8, 28.2, 27.9, 30.8]),
    ('pdb_4k3x_e_H18N10', [32.6, 33.4, 32.9, 34.0, 31.8, 28.0, 34.4, 28.0, 27.7, 31.2, 34.3]),
    ('pdb_4fqk_c_B',      [21.0, 22.2, 21.6, 21.6, 20.5, 21.0, 21.4, 21.2, 20.9, 21.9, 21.7, 21.1]),
    #                     4f3ze 1ru7a 1ruyh 3ztna 2wr7c 5hmge 2fk0k 4dj6a 1ti8a 1jsda 4i78b 4k3xe 4fqkc
])

# Z_SCORES will be a dict of dicts, indexed first by query id and then by
# subject id, holding the corresponding Z scores (taken from
# _Z_SCORES_TABLE).  We'll walk through _Z_SCORES_TABLE to fill in Z_SCORES.
Z_SCORES = {}

for index, queryId in enumerate(_Z_SCORES_TABLE):
    Z_SCORES[queryId] = {}
    for subjectId in _Z_SCORES_TABLE:
        if queryId != subjectId:
            try:
                score = _Z_SCORES_TABLE[subjectId][index]
            except IndexError:
                # The logic here may appear backwards, but it's not. To
                # convince yourself either think about it or look at the
                # testLogic test in test/performance/test_polymerase.py
                score = Z_SCORES[subjectId][queryId]

            Z_SCORES[queryId][subjectId] = score
