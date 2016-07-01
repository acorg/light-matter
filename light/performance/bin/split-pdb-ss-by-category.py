#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict

from dark.fasta_ss import SSFastaReads
from dark.reads import SSAAReadWithX, Reads

import sys
import argparse

parser = argparse.ArgumentParser(
    description=(
        'Read PDB secondary structure from stdin and split it into separate '
        'FASTA files according to specified category information. The '
        'sequence ids found on stdin be in the form e.g., pdb_2hla_a as '
        'produced by clean-pdb-ss-fasta.py.'))

parser.add_argument(
    '--prefix', default='pdb-category-',
    help=('A string filename prefix to use for output files. One output '
          'file will be created for each category. Output file names will '
          'be the prefix followed by the category name with a .fasta '
          'suffix.'))

parser.add_argument(
    '--categories', required=True,
    help=('A file specifying which PDB sequences go into which category. The '
          'file should contain space separated fields. The first field is '
          'the category name and subsequent fields (if any) are PDB sequence '
          'ids. Sequence ids are mapped to lower case.'))

parser.add_argument(
    '--keepChain', default=False, action='store_true',
    help=('If True, the chain information on PDB sequence ids will be '
          'considered in matching. If not, the chain is ignored, making '
          'it simple to extract all chains via a single sequence id.'))

parser.add_argument(
    '--allowSeqsInMultipleCategories', default=False, action='store_true',
    help=('If True, sequence ids may be in multiple categories. If not, '
          'a ValueError is raised if a sequence id appears in more than '
          'one category.'))

parser.add_argument(
    '--ignoreUncategorizedSequences', default=False, action='store_true',
    help=('If True, sequence ids that are not in any category will be '
          'reported to stderr. If not, they are silently ignored.'))

args = parser.parse_args()

# sequenceIdToCategories is keyed by a sequence id and its values hold all
# the categories that sequence id is in.
sequenceIdToCategories = defaultdict(set)

# Hold the name of all categories seen in args.categories.
categoryNames = set()

# Read all categories and their sequence ids.

with open(args.categories) as fp:
    for line in fp:
        fields = line.split()
        category = fields.pop(0)
        sequenceIds = fields
        if category in categoryNames:
            raise ValueError('Duplicate category %r found in %r' %
                             (category, args.categories))
        categoryNames.add(category)
        for sequenceId in map(str.lower, sequenceIds):
            if sequenceId in sequenceIdToCategories:
                if category in sequenceIdToCategories[sequenceId]:
                    # The sequence id has been given more than once in the
                    # *same* category.  We could just ignore this, but this
                    # kind of thing is often indicative of an error in
                    # whatever was used to produce the input, so I prefer
                    # to flag it.
                    raise ValueError(
                        'Sequence id %r found more than once in category %r '
                        'in %r' % (sequenceId, category, args.categories))
                else:
                    if not args.allowSeqsInMultipleCategories:
                        sequenceIdToCategories[sequenceId].add(category)
                        raise ValueError(
                            'Sequence id %r found in multiple categories (%s) '
                            'in %r' %
                            (sequenceId,
                             ', '.join(sorted(
                                 sequenceIdToCategories[sequenceId])),
                             args.categories))
            sequenceIdToCategories[sequenceId].add(category)

# Read the PDB sequence information and add each sequence to its category.
#
# Sequence ids must be in the form e.g., pdb_2hla_a as produced by the
# clean-pdb-ss-fasta.py (in this directory).

sequencesByCategory = defaultdict(Reads)

for sequence in SSFastaReads(sys.stdin, readClass=SSAAReadWithX,
                             checkAlphabet=0):
    pdb, sequenceId, chain = sequence.id.split('_')
    assert pdb == 'pdb' and len(chain) == 1, (
        'Unrecognized PDB id %r found on stdin.' % sequence.id)
    if args.keepChain:
        sequenceId += '_' + chain

    if sequenceId in sequenceIdToCategories:
        for category in sequenceIdToCategories[sequenceId]:
            sequencesByCategory[category].add(sequence)
    else:
        if not args.ignoreUncategorizedSequences:
            print('Sequence %r on stdin is not in any category.' % sequence.id,
                  file=sys.stderr)

# Write out all the sequences for each category. Loop over the keys of
# categoryNames (not sequencesByCategory) to make sure we write an empty
# file for categories that have no sequences in them.

for category in sorted(categoryNames):
    filename = args.prefix + category + '.fasta'
    sequences = sequencesByCategory[category]
    with open(filename, 'w') as fp:
        for sequence in sequences:
            print(sequence.toString(), file=fp, end='')
    print('Wrote %d sequences to %r' % (len(sequences), filename),
          file=sys.stderr)
