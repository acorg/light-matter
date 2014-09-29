#!/usr/bin/env python

import sys

from light.statistics import runStatistic, find

from dark.fasta import FastaReads

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s fastaFilename' % sys.argv[0]

    else:
        fastaFile = sys.argv[1]
        fastaReads = FastaReads(fastaFile)
        statistics = find()
        for statistic in statistics:
            count = runStatistic(statistic, fastaReads)
            print '%s: %d' % (statistic.NAME, count)
