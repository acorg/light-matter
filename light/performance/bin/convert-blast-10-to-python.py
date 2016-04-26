#!/usr/bin/env python

from __future__ import print_function

"""
Read BLAST format 10 (see below for details) input on stdin and write JSON
to stdout.

The input data MUST be the result of running:

    $ blastp -query queries.fasta -outfmt '10 qseqid sseqid bitscore' ...

That output format creates CSV with 3 fields per line.

Usage:

    $ convert-blast-10-to-json.py --subjects subjects.fasta \
        [--queries queries.fasta] < BLAST-output

If -queries is omitted, it will be assumed that the BLAST database was queried
with its subject sequences (i.e., that queries and subjects are identical).

Note that there doesn't seem to be a way to tell BLAST to only output one
score per query/subject match (-max_hsps limits the number of times a subject
is hit overall, not the number of times a query hits a subject). For this
reason the BLAST output may contain multiple lines for a given query/subject
combination. The code below uses the maximum bit score in this case.
"""

import sys
import csv
import argparse
from json import dumps

from dark.reads import AAReadWithX
from dark.fasta import FastaReads


parser = argparse.ArgumentParser(description='Find features in FASTA.')

parser.add_argument(
    '--subjects', required=True,
    help=('Specify the FASTA file used as the subjects to make the BLAST '
          'database to create the input.'))

parser.add_argument(
    '--queries',
    help=('Specify the FASTA file used as the queries to BLAST to create the '
          'input.'))

parser.add_argument(
    '--allowAsymmetricScores', default=False, action='store_true',
    help=('If True, the score of sequence A when used as a query against '
          'sequence B will not be coerced to match the score when B is '
          'matched against A. If False, both scores will be set to the '
          'maximum of the two match values.'))

parser.add_argument(
    '--compact', default=False, action='store_true',
    help=('If given, the output JSON will be compact (i.e., much harder for a '
          'human to read).'))

args = parser.parse_args()

subjects = dict((read.id, read) for read in
                FastaReads(args.subjects, readClass=AAReadWithX))

if args.queries:
    queries = dict((read.id, read) for read in
                   FastaReads(args.queries, readClass=AAReadWithX))
else:
    queries = subjects


result = {}

# Set all scores to zero.
for queryId in queries:
    result[queryId] = dict.fromkeys(subjects, 0.0)

for queryId, subjectId, value in csv.reader(sys.stdin):
    if queryId not in queries:
        raise RuntimeError(
            'Query %s found on stdin but not in %s' %
            (queryId, args.queries or args.subjects), file=sys.stderr)
    if subjectId not in subjects:
        raise RuntimeError(
            'Subject %r found on stdin but not in %r' %
            (subjectId, args.subjects), file=sys.stderr)
    value = float(value.strip())
    # Store the score if it's higher than what we already have.
    if value > result[queryId][subjectId]:
        result[queryId][subjectId] = value

if not args.allowAsymmetricScores:
    # Ensure query/subject and subject/query scores are equal. Compare the
    # reads to make sure their sequences (not just their read ids) are the
    # same.
    for queryId in queries:
        for subjectId in subjects:
            if (subjectId in queries and queryId in subjects and
                    queries[subjectId] == subjects[subjectId] and
                    queries[queryId] == subjects[queryId]):
                result[queryId][subjectId] = result[subjectId][queryId] = max(
                    result[queryId][subjectId], result[subjectId][queryId])


def adjustId(readId):
    """
    Change a read's id to use - instead of :
    """
    if readId.endswith(':sequence'):
        readId = readId[:-9]
    return readId.replace(':', '-')

new = {}
for queryId in queries:
    newQueryId = adjustId(queryId)
    assert newQueryId not in new
    new[newQueryId] = {}
    for subjectId in subjects:
        newSubjectId = adjustId(subjectId)
        assert newSubjectId not in new[newQueryId]
        new[newQueryId][newSubjectId] = result[queryId][subjectId]

print("""\
# THIS FILE IS AUTOMATICALLY GENERATED - DO NOT EDIT!
# See the Makefile for creation details.

BIT_SCORES = """, end='')

if args.compact:
    print(dumps(new, separators=(',', ':')))
else:
    print(dumps(new, indent='  ', sort_keys=True))
