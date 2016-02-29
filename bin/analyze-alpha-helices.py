#!/usr/bin/env python

"""
Read FASTA from stdin and examine it for features (landmarks and trig points).
Write a summary of the number of features that are found in the input for each
of the feature finders used.

The idea here is that you can provide input that contains known features (e.g.,
taken from PDB secondary structures) and check to see if our finders can find
the features.
"""

import sys
import argparse

from dark.aa import PROPERTY_CLUSTERS
from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.reads import AARead, SSAARead

from light.landmarks import THAlphaHelix
from light.performance.stats import Stats

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

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help=('If True, print detailed information.'))


def main(args):
    verbose = args.verbose
    minSequenceLength = args.minSequenceLength

    PEAK_HYDROPATHY = THAlphaHelix.PEAK_HYDROPATHY
    TROUGH_HYDROPATHY = THAlphaHelix.TROUGH_HYDROPATHY

    if args.pdb:
        readClass = SSAARead
        readsClass = SSFastaReads
    else:
        readClass = AARead
        readsClass = FastaReads

    readCount = 0

    length = Stats('Alpha helix length')
    initial = Stats('Initial non-peak non-trough AAs')
    final = Stats('Final non-peak non-trough AAs')
    extrema = Stats('Number of extrema')
    consecutivePeaks = Stats('Consecutive peaks')
    consecutiveTroughs = Stats('Consecutive troughs')
    noExtremaLength = Stats('Length of helices with no extrema')

    for read in readsClass(sys.stdin, readClass=readClass, checkAlphabet=0):
        readLen = len(read)

        if readLen < minSequenceLength:
            continue

        readCount += 1
        length.add(readLen)

        if verbose:
            print('read %d: len=%d %s id=%s seq=%r' % (
                readCount, readLen, read.id, read.sequence))

        # States.
        LOOKING = 0
        AWAITING_PEAK = 1
        AWAITING_TROUGH = 2

        state = LOOKING

        for offset, aa in enumerate(read.sequence):
            try:
                hydropathy = PROPERTY_CLUSTERS[aa]['hydropathy']
            except KeyError:
                hydropathy = 0

            if state == LOOKING:
                if hydropathy >= PEAK_HYDROPATHY:
                    state = AWAITING_TROUGH
                    finalExtremaOffset = offset
                    extremaCount = 1
                    consecutivePeakCount = 1
                    consecutiveTroughCount = 0
                    initial.add(offset)
                elif hydropathy <= TROUGH_HYDROPATHY:
                    state = AWAITING_PEAK
                    finalExtremaOffset = offset
                    extremaCount = 1
                    consecutivePeakCount = 0
                    consecutiveTroughCount = 1
                    initial.add(offset)

            elif state == AWAITING_PEAK:
                if hydropathy >= PEAK_HYDROPATHY:
                    finalExtremaOffset = offset
                    state = AWAITING_TROUGH
                    extremaCount += 1
                    consecutivePeakCount = 1
                    consecutiveTroughs.add(consecutiveTroughCount)
                    consecutiveTroughCount = 0
                elif hydropathy <= TROUGH_HYDROPATHY:
                    finalExtremaOffset = offset
                    consecutiveTroughCount += 1
                    extremaCount += 1

            elif state == AWAITING_TROUGH:
                if hydropathy <= TROUGH_HYDROPATHY:
                    finalExtremaOffset = offset
                    state = AWAITING_PEAK
                    extremaCount += 1
                    consecutivePeaks.add(consecutivePeakCount)
                    consecutivePeakCount = 0
                    consecutiveTroughCount = 1
                elif hydropathy >= PEAK_HYDROPATHY:
                    finalExtremaOffset = offset
                    consecutivePeakCount += 1
                    extremaCount += 1

        if state == LOOKING:
            noExtremaLength.add(readLen)
        else:
            final.add(readLen - finalExtremaOffset)
            extrema.add(extremaCount)

    print('Processed %d alpha helices.' % readCount)

    print(length)
    print(extrema)
    print(initial)
    print(final)
    print(consecutivePeaks)
    print(consecutiveTroughs)
    print(noExtremaLength)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
