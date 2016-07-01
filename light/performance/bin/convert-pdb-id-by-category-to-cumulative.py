#!/usr/bin/env python

"""
Read (from stdin) information giving categories of sequences that entered PDB.
Write (to stdout) the same information but in a cumulative way. E.g., given:

1996 1aaa 2bbb 3ccc
1997 4ddd 3eee 2fff 4ggg
1998 3hhh 1iii

We would print:

1996-1996 1aaa 2bbb 3ccc
1996-1997 1aaa 2bbb 3ccc 4ddd 3eee 2fff 4ggg
1996-1998 1aaa 2bbb 3ccc 4ddd 3eee 2fff 4ggg 3hhh 1iii

We use this script with ../../../data/pdb-structures-by-year.txt (per-year
PDB sequences) on input and use the output with split-pdb-ss-by-category.py
to break up the PDB ss.txt file in a cumulative way.

Note that although this script (and split-pdb-ss-by-category.py) are written in
terms of categories, in practice (at the moment) this means "years". Because
the code doesn't insist that categories could be years or numeric etc., you
need to give the per-category input in the order you want it accumulated. In
the case of years this would normally be going forward in time. Other
categories and accumulation orders can be imagined.
"""

from __future__ import print_function

import sys
import argparse

parser = argparse.ArgumentParser(
    description=('Read (from stdin) information giving categories of '
                 'sequences in PDB. Write (to stdout) the same information '
                 'but in a cumulative way. See the source code for an '
                 'example.'))

parser.add_argument(
    '--omitSummary', default=False, action='store_true',
    help=('If specified, the summary of how many sequence ids are in each '
          'cumulative category will not be printed.'))

args = parser.parse_args()

cumulative = set()
categories = set()
first = True

for line in sys.stdin:
    fields = line.split()
    category = fields.pop(0)
    if category in categories:
        raise ValueError('Duplicate category %r found on stdin' % category)
    if first:
        firstCategory = category
        first = False
    sequenceIds = set(fields)
    if cumulative & sequenceIds:
        raise ValueError(
            'Category %r contains sequences (%s) that had already appeared on '
            'stdin' % (category, ', '.join(sorted(cumulative & sequenceIds))))
    cumulative |= sequenceIds
    cumulativeCategory = firstCategory + '-' + category
    if not args.omitSummary:
        print('Category %s has %d sequences.' %
              (cumulativeCategory, len(cumulative)), file=sys.stderr)

    print('%s %s' % (cumulativeCategory, ' '.join(sorted(cumulative))))
