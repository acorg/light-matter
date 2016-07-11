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
from light.performance.pdb import loadObsolete, loadResolution

_DATA_DIR = dirname(dirname(light.__file__))
_DEFAULT_PDB_DELETION_FILE = join(
    _DATA_DIR, 'data', 'pdb-20160711-obsolete.txt')
_DEFAULT_PDB_RESOLUTION_FILE = join(
    _DATA_DIR, 'data', 'pdb-20160711-resolution.txt')

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
    '--pdbResolutionFile', default=_DEFAULT_PDB_RESOLUTION_FILE,
    metavar='PDB-resolution-file',
    help=('A file containing the resolution PDB structures, such as that '
          'found at ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx '
          'Specify /dev/null to skip filtering on PDB resolution.'))

parser.add_argument(
    '--maxResolution', default=3.0, type=float,
    help=('The maximum resolution (in Angstroms) PDB sequence to keep. Only '
          'sequences with this resolution or less will be kept.'))

parser.add_argument(
    '--minLength', default=4, type=int,
    help='The minimum length sequence to keep. Use 0 to keep all sequences.')

parser.add_argument(
    '--printIgnoredIds', action='store_true', default=False,
    help=('If given, print the PDB ids of sequences on input that are '
          'ignored due to being too short or their presence in the '
          '--pdbDeletionFile file'))

parser.add_argument(
    '--discardNMR', action='store_true', default=False,
    help=('If given, structures determined by NMR will be discarded.'))


args = parser.parse_args()

pdbDeletions = set(loadObsolete(args.pdbDeletionFile))
pdbResolutions = loadResolution(args.pdbResolutionFile)

deleted = set()
tooShort = set()
poorResolution = set()
nmr = set()

minLength = args.minLength
maxResolution = args.maxResolution
discardNMR = args.discardNMR

# The NMR resolution is assigned to PDB structures that were obtained via
# NMR (as opposed to crystallization). For now we keep all such
# structures. This is briefly mentioned at
# http://www.rcsb.org/pdb/static.do?p=general_information/about_pdb/\
# summaries.html
NMR_RESOLUTION = -1.0

for sequence in SSFastaReads(sys.stdin, checkAlphabet=0):
    pdbId = sequence.id.split(':', maxsplit=1)[0]

    # Check if this sequence has been deleted from PDB.
    if pdbDeletions:
        if pdbId in pdbDeletions:
            deleted.add(pdbId)
            continue

    # Check if the resolution on this sequence is good (i.e., numerically
    # low) enough.
    if pdbResolutions:
        try:
            resolution = pdbResolutions[pdbId]
        except KeyError:
            print('PDB id %r has unknown resolution!' % pdbId, file=sys.stderr)
        else:
            if resolution == NMR_RESOLUTION and discardNMR:
                nmr.add(pdbId)
                continue
            if resolution > maxResolution:
                poorResolution.add(pdbId)
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
        print('%d sequence ids ignored due to being marked as obsolete in '
              'PDB: %s.' % (len(deleted), ', '.join(sorted(deleted))),
              file=sys.stderr)
    else:
        print('%d sequence ids ignored due to being PDB marked as obsolete in '
              'PDB.' % len(deleted), file=sys.stderr)

if tooShort:
    if args.printIgnoredIds:
        print('%d sequences were ignored due to being too short (< %d '
              'residues): %s.' % (
                  len(tooShort), minLength, ', '.join(sorted(tooShort))),
              file=sys.stderr)
    else:
        print('%d sequences were ignored due to being too short (< %d '
              'residues).' % (len(tooShort), minLength), file=sys.stderr)

if poorResolution:
    if args.printIgnoredIds:
        print('%d sequences were ignored due to poor resolution (> %.3f '
              'angstroms): %s.' % (len(poorResolution), maxResolution,
                                   ', '.join(sorted(poorResolution))),
              file=sys.stderr)
    else:
        print('%d sequences were ignored due to poor resolution (> %.3f '
              'angstroms).' %
              (len(poorResolution), maxResolution), file=sys.stderr)

if nmr:
    if args.printIgnoredIds:
        print('%d sequences were ignored as they were determined by NMR: %s.'
              % (len(nmr), ', '.join(sorted(nmr))),
              file=sys.stderr)
    else:
        print('%d sequences were ignored as they were determined by NMR. ' %
              len(nmr), file=sys.stderr)

if (deleted or tooShort or poorResolution or nmr) and not args.printIgnoredIds:
    print('Use --printIgnoredIds to print ignored seqeunce ids.',
          file=sys.stderr)
