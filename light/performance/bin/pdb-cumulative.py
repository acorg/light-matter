#!/usr/bin/env python

"""
"""

from __future__ import print_function, division

import time
import os
import re
from os.path import join
import argparse

from dark.fasta_ss import SSFastaReads

from light.landmarks import findLandmark
from light.performance.evaluate import evaluateMatchNoPrefix


def getSubstrings(length, sequence):
    """
    Take a sequence and make substrings of a specified length.

    @param length: An C{int} length of the desired subsequences.
    @param sequence: A C{str} amino acid sequence.
    @return: A generator that yields C{str} substrins of C{sequence}
        type as C{sequence}).
    """
    seqLen = len(sequence)
    if length == seqLen:
        yield sequence
    elif length < seqLen:
        numberOfSubsequences = seqLen - length + 1
        for i in range(numberOfSubsequences):
            yield sequence[i:length + i]


def tpfp(pdb, substring, structureType):
    """
    """
    truePositive = falsePositive = 0
    regex = re.compile(substring)
    for read in pdb:
        structure = read.structure
        for match in regex.finditer(read.sequence):
            start, end = match.start(), match.end()
            if evaluateMatchNoPrefix(structure, start, end,
                                     structureType):
                truePositive += 1
            else:
                falsePositive += 1

    return [truePositive, falsePositive]


def tpfpKey(substring):
    """
    """
    tp, fp = tpfpUpToPreviousYear[substring]
    return tp / (tp + fp), tp, -fp

startTime = time.time()

parser = argparse.ArgumentParser(
    description=('Cumulative analysis of PDB'))

parser.add_argument(
    '--featureType', default='PDB AlphaHelix',
    choices=('PDB AlphaHelix', 'PDB AlphaHelix_3_10',
             'PDB AlphaHelix_combined', 'PDB AlphaHelix_pi',
             'PDB ExtendedStrand'),
    help='The type of structure to analyze.')

parser.add_argument(
    '--minLength', type=int, default=4,
    help='The minimum length subsequence to produce.')

parser.add_argument(
    '--maxLength', type=int,
    help=('The maximum length subsequence to produce. If unspecified, '
          'all subsequences at least as long as minLength will be produced.'))

parser.add_argument(
    '--sort', action='store_true', default=False,
    help='If true, sort output files by tp/(tp+fp) score and by substring.')

args = parser.parse_args()

featureType = args.featureType
assert(featureType.startswith('PDB '))
finder = findLandmark(featureType)()
ahType = featureType[4:]
structureLetter = finder.STRUCTURE_LETTER

# Years to operate on (inclusive).
firstYear = 1972
lastYear = 1995

years = list(map(str, range(firstYear, lastYear + 1)))

pdbUpToPreviousYear = []
tpfpUpToPreviousYear = {}

for year in years:
    yearStartTime = time.time()
    print('Processing year %s' % year)
    fastaFile = join(year, 'pdb-' + year + '.fasta')

    pdbThisYear = list(SSFastaReads(fastaFile, checkAlphabet=0))
    tpfpThisYear = {}

    print('  Read %d sequence%s from %s' % (
        len(pdbThisYear), '' if len(pdbThisYear) == 1 else 's', fastaFile))

    for read in pdbThisYear:
        for feature in finder.find(read):
            start = feature.offset
            end = start + feature.length
            sequence = read.sequence[start:end]
            if args.maxLength is None:
                maxLength = len(sequence)
            else:
                maxLength = min(args.maxLength, len(sequence))
            for length in range(args.minLength, maxLength + 1):
                for substring in getSubstrings(length, sequence):
                    # If not already calculated, work out the true and
                    # false positive counts for this substring.

                    # Previous years:
                    if substring not in tpfpUpToPreviousYear:
                        tpfpUpToPreviousYear[substring] = tpfp(
                            pdbUpToPreviousYear, substring, structureLetter)

                    # This year:
                    if substring not in tpfpThisYear:
                        tpfpThisYear[substring] = tpfp(
                            pdbThisYear, substring, structureLetter)

    print('  %d new %s substring%s occurred in %s' % (
        len(tpfpThisYear), ahType, '' if len(tpfpThisYear) == 1 else 's',
        year))

    # For all substrings that have been seen in previous years but which
    # are not in this year, evaluate them against this year's
    # sequences. This will cause the false positive counts for some of them
    # to increase.
    for substring in tpfpUpToPreviousYear:
        if substring not in tpfpThisYear:
            tp, fp = tpfp(pdbThisYear, substring, structureLetter)
            assert tp == 0, (
                "Oops! Substring %r should not have any true positives!" %
                substring)
            prev = tpfpUpToPreviousYear[substring]
            prev[0] += tp
            prev[1] += fp

    # Merge this year's tp/fp counts into those for all previous years.
    for substring, (tp, fp) in tpfpThisYear.items():
        prev = tpfpUpToPreviousYear[substring]
        prev[0] += tp
        prev[1] += fp

    print('  %d substring%s have been seen through %s' % (
        len(tpfpUpToPreviousYear),
        '' if len(tpfpUpToPreviousYear) == 1 else 's', year))

    # Add this year's PDB to the previous PDB.
    pdbUpToPreviousYear.extend(pdbThisYear)

    print('  PDB through %s contains %d sequence%s' % (
        year, len(pdbUpToPreviousYear),
        '' if len(pdbUpToPreviousYear) == 1 else 's'))

    # Write out the tp/fp results to date.
    resultDir = join(year, ahType)
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)
    resultFile = join(resultDir, year + '.out')

    if args.sort:
        # Sort first by the substring (ascending), then by tp/(tp+fp) score
        # (reversed).
        substrings = sorted(tpfpUpToPreviousYear)
        substrings.sort(key=tpfpKey, reverse=True)
    else:
        substrings = tpfpUpToPreviousYear

    with open(resultFile, 'w') as resultFp:
        for substring in substrings:
            tp, fp = tpfpUpToPreviousYear[substring]
            print('%s %d %d' % (substring, tp, fp), file=resultFp)

    print('  %s done in %.2f seconds' % (year, time.time() - yearStartTime))

print('Overall run took %.2f seconds' % (time.time() - startTime))
