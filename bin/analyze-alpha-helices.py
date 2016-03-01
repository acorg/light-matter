#!/usr/bin/env python

"""
Read FASTA containing known alpha helices from stdin (likely obtained
by running extract-alpha-helice-from-pdb-ss.py) and compute statistics on
their properties (primarily amino acid hydropathy). Write a summary to
stdout.
"""

import sys
import argparse

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.reads import AARead, SSAARead

from light.performance.alpha_helix import analyzeAlphaHelices

parser = argparse.ArgumentParser(description='Find features in FASTA.')

parser.add_argument(
    '--minSequenceLength', '--msl', type=int, default=0,
    help=('Sequences in the input whose length is less than this value '
          'will not be considered.'))

parser.add_argument(
    '--pdb', default=False, action='store_true',
    help=('If True, the input is in PDB format, in which two FASTA records '
          'are used per sequence. The first holds the AA sequence and the '
          'second has the structure. If False, the input is taken as regular '
          'FASTA.'))


def main(args):
    """
    Use analyzeAlphaHelices to get a collection of alpha helix reads analyzed
    and print the results.

    @param args: Command line arguments as returned by the C{argparse}
        C{parse_args} method.
    """
    if args.pdb:
        readClass = SSAARead
        readsClass = SSFastaReads
    else:
        readClass = AARead
        readsClass = FastaReads

    reads = readsClass(sys.stdin, readClass=readClass, checkAlphabet=0)

    if args.minSequenceLength:
        reads = reads.filter(minLength=args.minSequenceLength)

    analysis = analyzeAlphaHelices(reads)

    print('Processed %d alpha helices.' % len(reads))

    for key in ('length', 'extremaCount', 'initial', 'final',
                'consecutivePeaks', 'consecutiveTroughs', 'noExtremaLength'):
        if analysis[key]:
            print(analysis[key])
        else:
            print('No values recorded for', key)


if __name__ == '__main__':
    main(parser.parse_args())
