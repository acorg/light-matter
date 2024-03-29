#!/usr/bin/env python

from __future__ import print_function

import argparse
import re
import sys

from dark.fasta_ss import SSFastaReads
from dark.fasta import FastaReads
from dark.reads import AAReadWithX

from light.performance.evaluate import evaluateMatch, evaluateMatchNoPrefix


parser = argparse.ArgumentParser(
    description=('For each helix in the file given on stdin, evaluate its '
                 'true and false positives and print the result to stdout.'))

parser.add_argument(
    '--pdbFile', help='A filename of the pdb file containing sequences and '
    'their structure annotation.')

parser.add_argument(
    '--evaluateNoPrefix', default=True, type=bool,
    help=('If True the evaluateMatchNoPrefix function will be used to '
          'evaluate the helix. If False use the evaluateMatch function.'))

parser.add_argument(
    '--structureType', default='H', choices={'H', 'G', 'I', 'E', 'K'},
    help=('The type of structure that should be evaluated against. '
          'H: Alpha helix, G: Alpha helix 3 10, I: Alpha helix pi, E: '
          'Extended strand, K: Combined alpha helix.'))

args = parser.parse_args()

pdbReads = [(read.sequence, read.structure) for read in
            SSFastaReads(args.pdbFile, checkAlphabet=0)]
helices = FastaReads(sys.stdin, readClass=AAReadWithX, checkAlphabet=0)

if args.evaluateNoPrefix:
    evaluationFunction = evaluateMatchNoPrefix
else:
    evaluationFunction = evaluateMatch

for i, helix in enumerate(helices):
    truePositive = falsePositive = 0
    helixSequence = helix.sequence
    if ('X' in helixSequence or 'Z' in helixSequence or 'B' in
            helixSequence):
        continue
    else:
        uniqueRegex = re.compile(helixSequence)
        for sequence, structure in pdbReads:
            for match in uniqueRegex.finditer(sequence):
                start = match.start()
                end = match.end()
                if evaluationFunction(structure, start, end,
                                      args.structureType):
                    truePositive += 1
                else:
                    falsePositive += 1

        print('%s %d %d' % (helix.sequence, truePositive, falsePositive))
