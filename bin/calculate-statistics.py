#!/usr/bin/env python

import sys

from light.statistics import runStatistic
from light.statistics.all_bases import HasAllBases

from dark.fasta import FastaReads

statistics = Statistic.__subclasses__()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s fastaFilename' % sys.argv[0]

    else:
        fastaFile = sys.argv[1]
        fastaReads = FastaReads(fastaFile)
        for statistic in statistics:
            count = runStatistic(statistic, fastaReads)
            print '%s: %d' % (statistic.NAME, count)
