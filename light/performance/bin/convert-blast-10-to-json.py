#!/usr/bin/env python

from __future__ import print_function

"""
Read BLAST format 10 (see below for details) input on stdin and write JSON
to stdout.

The input data MUST be the result of running:

    $ blastp -query queries.fasta -outfmt '10 qseqid sseqid bitscore' ...

That output format creates CSV with 3 fields per line.

If you created the BLAST database using the queries.fasta file, you will
probably find it useful to use --fasta queries.fasta on the command line.
This will ensure that the output contains a 0.0 score when a query doesn't
match a subject. Queries that don't match any subjects will also be in the
output when using --fasta.

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
    '--fasta',
    help=('Specify a FASTA file that was used as both query and subjects in '
          'the BLAST run that created the input. This will be read to get all '
          'sequence ids, to make it possible to output a JSON object that '
          'includes keys and values for query and subject ids that were not '
          'matched by anything.'))

parser.add_argument(
    '--allowAsymmetricScores', default=False, action='store_true',
    help=('If True, the score of sequence A when used as a query against '
          'sequence B will not be coerced to match the score when B is '
          'matched against A. If False, both scores will be set to the '
          'maximum of the two match values.'))

parser.add_argument(
    '--includeSelfMatches', default=False, action='store_true',
    help='If True, bit scores of a sequence against itself will be included.')

parser.add_argument(
    '--omitMissingScores', default=False, action='store_true',
    help=('If given, bit scores of 0.0 will not be included for queries that '
          'don\'t match subjects.'))

parser.add_argument(
    '--compact', default=False, action='store_true',
    help=('If given, the output JSON will be compact (i.e., much harder for a '
          'human to read).'))

args = parser.parse_args()

if args.fasta:
    allIds = set([read.id for read in
                  FastaReads(args.fasta, readClass=AAReadWithX)])
else:
    allIds = None

result = {}

subjectIds = set()

reader = csv.reader(sys.stdin)
for queryId, subjectId, value in reader:
    subjectIds.add(subjectId)
    if queryId != subjectId or args.includeSelfMatches:
        value = float(value.strip())
        if queryId in result:
            if subjectId in result[queryId]:
                # Store the new score if it's higher than what we already have.
                if value > result[queryId][subjectId]:
                    result[queryId][subjectId] = value
            else:
                result[queryId][subjectId] = value
        else:
            result[queryId] = {subjectId: value}

if not args.allowAsymmetricScores:
    # Ensure query/subject and subject/query scores both exist and are equal.
    for queryId in list(result):
        for subjectId in result[queryId]:
            if subjectId not in result:
                result[subjectId] = {}
            if queryId in result[subjectId]:
                score = max(result[queryId][subjectId],
                            result[subjectId][queryId])
                result[queryId][subjectId] = result[subjectId][queryId] = score
            else:
                result[subjectId][queryId] = result[queryId][subjectId]

if not args.omitMissingScores:
    # Add 0.0 scores for query/subject pairs that are not hit by reads.
    for queryId in allIds or result:
        if queryId not in result:
            result[queryId] = {}
        for subjectId in allIds or subjectIds:
            if queryId != subjectId or args.includeSelfMatches:
                if subjectId not in result[queryId]:
                    result[queryId][subjectId] = 0.0

if args.compact:
    print(dumps(result, separators=(',', ':')))
else:
    print(dumps(result, indent='  ', sort_keys=True))
