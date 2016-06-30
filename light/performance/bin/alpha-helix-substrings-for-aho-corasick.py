#!/usr/bin/env python

from __future__ import print_function, division


import sys
import argparse

from light.performance.alpha_helix import selectSubstringsForAhoCorasick


def printSummary(result):
    """
    Print a result / processing summary.

    @param result: A C{dict} as returned by C{selectSubstringsForAhoCorasick}.
    """
    inputCount = result['inputCount']
    print('Kept %d of %d (%.2f%%) candidate substring%s seen on input.' %
          (len(result['substrings']), inputCount,
           len(result['substrings']) / inputCount * 100.0,
           '' if inputCount == 1 else 's'), file=sys.stderr)

    notEnoughTruePositives = result['notEnoughTruePositives']
    print('%d substring%s did not meet the minimum true positive '
          'requirement (%d).' %
          (notEnoughTruePositives, '' if notEnoughTruePositives == 1 else 's',
           args.minTruePositives),
          file=sys.stderr)

    fractionTooLow = result['fractionTooLow']
    print('%d substring%s did not have a high enough true positive '
          'fraction (%f).' %
          (fractionTooLow, '' if fractionTooLow == 1 else 's',
           args.minTruePositiveFraction),
          file=sys.stderr)

    inferior = result['inferior']
    if inferior == 1:
        print('1 substring was inferior to (at least) one of its own '
              'substrings.', file=sys.stderr)
    else:
        print('%d substrings were inferior to (at least) one of their own '
              'substrings.' % inferior, file=sys.stderr)


parser = argparse.ArgumentParser(
    description=('Read (from stdin) PDB alpha helix substrings with their '
                 'true positive and false positive counts and print (to '
                 'stdout) a selected subset of the substrings for use in '
                 'the Aho Corasick alpha helix finder. Stdin should have '
                 'lines with 3 space-separated fields: an alpha helix '
                 'substring, its integer true positive count, its integer '
                 'false positive count.'))

parser.add_argument(
    '--minTruePositives', type=int, default=0,
    help=('The minimum number of true positives a substring must have to be '
          'considered.'))

parser.add_argument(
    '--minTruePositiveFraction', type=float, default=0.0,
    help=('The minimum ratio of true positives to overall (true + false) '
          'positives a substring must have to be considered.'))

parser.add_argument(
    '--maxSubstrings', type=int, default=-1,
    help=('The maximum number of substrings to return. By default no limit '
          'is applied.'))

parser.add_argument(
    '--printCounts', default=False, action='store_true',
    help=('If True, the true positive count, false positive count, and true '
          'positive count / (false positive count + true positive count) '
          'fraction will be printed after each substring.'))

parser.add_argument(
    '--printSummary', default=False, action='store_true',
    help=('If True, print a summary of substring processing to show how many '
          'substrings were considered and what their fates were.'))

args = parser.parse_args()

result = selectSubstringsForAhoCorasick(
    sys.stdin,
    args.minTruePositives,
    args.minTruePositiveFraction,
    args.maxSubstrings)

if args.printCounts:
    for substring, counts in result['substrings']:
        print('%s %d %d %f' % (substring, counts[0], counts[1], counts[2]))
else:
    for substring, _ in result['substrings']:
        print(substring)

if args.printSummary:
    printSummary(result)
