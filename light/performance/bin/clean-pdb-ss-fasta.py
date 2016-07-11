#!/usr/bin/env python

"""
Read the PDB ss.txt file on stdin and do the following:

  1) Change sequence id lines like ">101M:A:sequence" to ">pdb_101m_a"
  2) Change structure id lines like ">101M:A:secstr" to ">pdb_101m_a:structure"
  3) Omit any sequences whose ids are found in the PDB deletions file (if
     one is specified).
  4) Omit any sequences that are too short.

Write the modified output to stdout.

The input will typically be the ss.txt secondary structure file available at
http://www.rcsb.org/pdb/static.do?p=download/http/index.html  Note that that
file contains spaces in its structure lines and these are not dealt with
properly by from Bio import SeqIO so they must be converted to hyphens before
arriving on stdin to this script. So typical usage should be:

    $ gzcat ss.txt.gz | tr ' ' - | clean-pdb-ss-fasta.py
"""

from __future__ import print_function
import sys
import argparse
from os.path import dirname, join

from dark.fasta_ss import SSFastaReads

import light
from light.performance.pdb import loadObsoletePDB


_DEFAULT_PDB_DELETION_FILE = join(
    dirname(dirname(light.__file__)), 'data', 'pdb-20160711-obsolete.txt')

parser = argparse.ArgumentParser(
    description=('Clean a PDB ss.txt file. See source code for description '
                 'of cleaning steps.'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '--pdbDeletionFile', default=_DEFAULT_PDB_DELETION_FILE,
    metavar='PDB-deletion-file',
    help=('A file containing details of deleted PDB structures, such as that '
          'found at ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat '
          'Specify /dev/null to skip PDB deletion.'))

parser.add_argument(
    '--minLength', default=4, type=int,
    help='The minimum length sequence to keep. Use 0 to keep all sequences.')

parser.add_argument(
    '--printIgnoredIds', action='store_true', default=False,
    help=('If given, print the PDB ids of sequences on input that are '
          'ignored due to being too short or their presence in the '
          '--pdbDeletionFile file'))


args = parser.parse_args()

pdbDeletions = set(loadObsoletePDB(args.pdbDeletionFile))

deleted = set()
tooShort = set()

minLength = args.minLength

for sequence in SSFastaReads(sys.stdin, checkAlphabet=0):
    pdbId = sequence.id.split(':', maxsplit=1)[0]

    # Check if this sequence has been deleted from PDB.
    if pdbDeletions:
        if pdbId in pdbDeletions:
            deleted.add(pdbId)
            continue

    if minLength > 0 and len(sequence) < minLength:
        tooShort.add(pdbId)
        continue

    if sequence.id.endswith(':sequence'):
        sequence.id = sequence.id[:-9]
    else:
        raise RuntimeError('Unrecognized sequence id line %r' % sequence.id)

    print(sequence.toString(), end='')


if deleted:
    if args.printIgnoredIds:
        print('%d deleted sequence ids ignored: %s.' % (
            len(deleted), ', '.join(sorted(deleted))), file=sys.stderr)
    else:
        print('%d deleted sequence ids ignored. Use --printIgnoredIds to '
              'see them.' % len(deleted), file=sys.stderr)

if tooShort:
    if args.printIgnoredIds:
        print('%d sequences were ignored due to being too short: %s.' % (
            len(tooShort), ', '.join(sorted(tooShort))), file=sys.stderr)
    else:
        print('%d sequences were ignored due to being too short. '
              'Use --printIgnoredIds to see them.' %
              len(tooShort), file=sys.stderr)
