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

args = parser.parse_args()

sequenceIdToCategory = {}
categories = set()

# Read all categories and their sequence ids.

with open(args.categories) as fp:
    for line in fp:
        fields = line.split()
        category = fields.pop(0)
        sequenceIds = fields
        if category in categories:
            raise ValueError('Duplicate category %r found in %r' %
                             (category, args.categories))
        categories.add(category)
        for sequenceId in map(str.lower, sequenceIds):
            if sequenceId in sequenceIdToCategory:
                if sequenceIdToCategory[sequenceId] == category:
                    # At the moment we only allow a sequence to be in a
                    # single category. That could be relaxed.
                    raise ValueError(
                        'Sequence id %r found more than once in category %r '
                        'in %r' % (sequenceId, category, args.categories))
                else:
                    raise ValueError(
                        'Sequence id %r found in multiple categories (%r '
                        'and %r) in %r' %
                        (sequenceId, category,
                         sequenceIdToCategory[sequenceId],
                         args.categories))
            sequenceIdToCategory[sequenceId] = category

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

    try:
        category = sequenceIdToCategory[sequenceId]
    except KeyError:
        print('Sequence %r on stdin is not in any category.' % sequence.id,
              file=sys.stderr)
    else:
        sequencesByCategory[category].add(sequence)

# Write out all the sequences for each category. Loop over the keys of
# categories (not sequencesByCategory) to make sure we write an empty file
# for categories that have no sequences in them.

for category in sorted(categories):
    filename = args.prefix + category + '.fasta'
    sequences = sequencesByCategory[category]
    with open(filename, 'w') as fp:
        for sequence in sequences:
            print(sequence.toString(), file=fp, end='')
    print('Wrote %d sequences to %r' % (len(sequences), filename),
          file=sys.stderr)
