#!/usr/bin/env python

"""
Read the PDB ss.txt file on stdin and write (to stdout) lines that contain a
year followed by the PDB ids with accession dates for that year.

The input will typically be our cleaned up version of the ss.txt secondary
structure file available at
http://www.rcsb.org/pdb/static.do?p=download/http/index.html  A cleaned
version of this file is available in data/pdb-20160711-ss.txt.bz2
"""

from __future__ import print_function

import sys
import argparse
from os.path import dirname, join
from collections import defaultdict

from dark.fasta_ss import SSFastaReads

import light
from light.performance.pdb import loadEntries, pythonNameToPdbName

_DATA_DIR = dirname(dirname(light.__file__))
_DEFAULT_PDB_ENTRIES_FILE = join(
    _DATA_DIR, 'data', 'pdb-20160711-entries.txt')

parser = argparse.ArgumentParser(
    description=('Clean a PDB ss.txt file. See source code for description '
                 'of cleaning steps.'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '--pdbEntriesFile', default=_DEFAULT_PDB_ENTRIES_FILE,
    metavar='PDB-entries-file',
    help=('A file containing details of PDB structures, such as that '
          'found at ftp://ftp.wwpdb.org/pub/pdb/data/status/entries.dat'))

parser.add_argument(
    '--quiet', default=False, action='store_true',
    help=('If given, do not write a summary of how many sequence ids were '
          'processed and put into each category'))

args = parser.parse_args()

pdbEntries = loadEntries(args.pdbEntriesFile)

sequenceCount = yearFoundCount = 0
noYear = set()
yearToIds = defaultdict(set)

for sequence in SSFastaReads(sys.stdin, checkAlphabet=0):
    sequenceCount += 1
    pdbId = pythonNameToPdbName(sequence.id).split(':', maxsplit=1)[0]

    try:
        entry = pdbEntries[pdbId]
    except KeyError:
        noYear.add(pdbId)
    else:
        yearToIds[entry['year']].add(pdbId)
        yearFoundCount += 1

verbose = not args.quiet

noYearCount = len(noYear)

if verbose:
    print('%d sequence%s read, %d had a year, %d had no year.' %
          (sequenceCount, '' if sequenceCount == 1 else 's',
           yearFoundCount, noYearCount), file=sys.stderr)

for year in sorted(yearToIds):
    print('%d %s' % (year, ' '.join(sorted(yearToIds[year]))))
    if verbose:
        count = len(yearToIds[year])
        print('%d: %d sequence%s' % (year, count, '' if count == 1 else 's'),
              file=sys.stderr)

if noYear:
    print('%d sequence ids found on stdin were not present in %r: %s' %
          (noYearCount, args.pdbEntriesFile, sorted(noYear)),
          file=sys.stderr)
