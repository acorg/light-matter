#!/usr/bin/env python

"""
Read PDB secondary structure records from stdin, find their alpha
helices (from the predicted secondary structure) and print just the alpha
helices to stdout.

Optionally, drop the structure information from the output (to produce regular
FASTA) and/or add a margin of amino acids around the found alpha helices.

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

from light.landmarks import PDB_AlphaHelix


parser = argparse.ArgumentParser(
    description=('Extract alpha helix subsequences from predicted PDB '
                 'secondary structures.'))

parser.add_argument(
    '--margin', type=int, default=0,
    help=('The number of AA residues to keep on either side of the '
          'extracted helices. It is useful to keep some context in order '
          'to give other alpha helix finders (e.g., GOR4) that you might '
          'run on the extracted helices a chance to establish context. '
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

args = parser.parse_args()
finder = PDB_AlphaHelix()

dropStructure = args.dropStructure
margin = args.margin

if margin < 0:
    raise ValueError('Margin must be non-negative.')

# The ss.txt file available at PDB has sequences that contain (at least)
# 'X' and 'U'. So, for now, read it without checking sequence alphabet.

for read in SSFastaReads(sys.stdin, checkAlphabet=0):
    for helix in finder.findWithMargin(read, margin):

        # Drop the ':sequence' suffix from read ids and add information
        # about the (1-based) offsets at which this helix was found.
        start = helix.offset - margin
        end = helix.offset + helix.length + margin
        readId = read.id.replace(':sequence', '') + ':%d-%d' % (start + 1, end)

        if dropStructure:
            helixWithMargin = AARead(readId, read.sequence[start:end])
        else:
            helixWithMargin = read[start:end]
            helixWithMargin.id = readId

        print(helixWithMargin.toString(format_='fasta'), end='')
