#!/usr/bin/env python

"""
Read PDB secondary structure records from stdin, find the desired structure
(from the predicted secondary structure) and print just the desired structure
to stdout.

Optionally, drop the structure information from the output (to produce regular
FASTA) and/or add a margin of amino acids around the found structure.

See the docstring for dark.fasta_ss.SSFastaReads for information on the PDB
secondary structure file format. Note the IMPORTANT NOTE section there! You
will likely want to run the ss.txt file through "tr ' ' -" on the UNIX command
line to convert all spaces to hyphens, before providing it as input to this
program.
"""

import sys
import argparse

from dark.fasta_ss import SSFastaReads
from dark.reads import AARead

from light.landmarks import findLandmark


parser = argparse.ArgumentParser(
    description=('Extract structure subsequences from predicted PDB '
                 'secondary structures.'))

parser.add_argument(
    '--margin', type=int, default=0,
    help=('The number of AA residues to keep on either side of the '
          'extracted structures. It is useful to keep some context in order '
          'to give other structure finders (e.g., GOR4) that you might '
          'run on the extracted structures a chance to establish context. '
          'Note that if a feature does not have this many neighboring '
          'residues (due to being too close to the start or end of a '
          'sequence, it will not be output.'))

parser.add_argument(
    '--dropStructure', '--ds', default=False, action='store_true',
    help=('If True, only write the sequence to the output, and do not include '
          'the matching structure information. This makes the output be '
          'regular one-record-per-sequence FASTA. If False (the default), '
          'the output will match the (PDB) input format, with sequence '
          'and matching structure information written as successive FASTA '
          'records.'))

parser.add_argument(
    '--feature', default='PDB AlphaHelix',
    choices=('PDB AlphaHelix', 'PDB AlphaHelix_3_10', 'PDB AlphaHelix_pi',
             'PDB ExtendedStrand'),
    help='The type of structure that should be extracted.')

args = parser.parse_args()

finder = findLandmark(args.structureType)

dropStructure = args.dropStructure
margin = args.margin

if margin < 0:
    raise ValueError('Margin must be non-negative.')

# The ss.txt file available at PDB has sequences that contain (at least)
# 'X' and 'U'. So, for now, read it without checking sequence alphabet.

for read in SSFastaReads(sys.stdin, checkAlphabet=0):
    for feature in finder.findWithMargin(read, margin):

        # Drop the ':sequence' suffix from read ids and add information
        # about the (1-based) offsets at which this feature was found.
        start = feature.offset - margin
        end = feature.offset + feature.length + margin
        readId = read.id.replace(':sequence', '') + ':%d-%d' % (start + 1, end)

        if dropStructure:
            featureWithMargin = AARead(readId, read.sequence[start:end])
        else:
            featureWithMargin = read[start:end]
            featureWithMargin.id = readId

        print(featureWithMargin.toString(format_='fasta'), end='')
