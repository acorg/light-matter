#!/usr/bin/env python

import sys

from light.utils import HasAllBases
from dark.fasta import FastaReads

statistics = [HasAllBases()]

def runStatistic(statistic, fastaReads):
    count = 0
    for read in fastaReads:
        result = statistic.evaluate(read)
        if result:
            count += 1
    print statistic.NAME, count 


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, ('Usage: %s fastaFilename' % sys.argv[0])

    else:
        fastaFile = sys.argv[1]
        fastaReads = FastaReads(fastaFile)
        for statistic in statistics:
            runStatistic(statistic, fastaReads)


