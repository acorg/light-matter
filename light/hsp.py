from __future__ import print_function

import sys
from operator import attrgetter

from light.bin_score import histogramBinFeatures


def _printBin(bin_):
    """
    Print detail of the contents of a histogram bin.

    @param bin_: A C{light.histogram.Histogram} bin, which is a C{list} of
        C{dict}s.
    """
    print('Bin has %d items:' % len(bin_), file=sys.stderr)
    for i, hashInfo in enumerate(bin_, start=1):
        print('  Item %d:' % i, file=sys.stderr)
        for key, value in hashInfo.items():
            # The 16 below is the length of the longest key (subjectTrigPoint).
            print('    %16s: %s' % (key, value), file=sys.stderr)


def normalizeBin(bin_, queryLen):
    """Examine a bin and return information about where the query and the subject
    begin and end. Return a dict with keys that allow the query and the subject
    to be displayed relative to the subject orientation (i.e., with start <
    stop for both the read and the match). The returned query indices are
    offsets into the subject. I.e., they indicate where on the subject the
    query lies.

    The diagram below shows one possible way in which the query and the
    subject can align and the terminology used. Note that there are many
    other possibilities. 'S' is used to indicate a subject feature offset,
    and 'Q' indicates a query feature offset.

                           queryStartInSubject             queryEndInSubject
                           |                               |
                           |                               |
                           |   subjectBinStart   subjectBinEnd
                           |   |                 |         |
                           |   |                 |         |
    Subject:  .................SSS..SSSS.....SSSS.................
    Query:                 ....QQQ..QQQQ.....QQQQ..........
                               |                 |
                               |                 |
                               queryBinStart     queryBinEnd

    NOTE: the returned queryStartInSubject value may be negative.  We
    consider the subject sequence to start at offset 0.  So if the query
    string has sufficient additional residues before the start of the
    alignment match, it may protrude to the left of the subject. Similarly,
    the returned queryEndInSubject can be greater than the subjectEnd.

    Here are two alignment examples that came up in the debugging of
    https://github.com/acorg/light-matter/issues/493
    Both of these can only happen (as far as we know!) when distanceBase > 1.0.

    The first illustrates how queryEndInSubject can be less than subjectBinEnd:

                           queryStartInSubject  queryEndInSubject
                           |                    |
                           |                    |
                           |   subjectBinStart  |subjectBinEnd
                           |   |                ||
                           |   |                ||
    Subject:  .................SSS..SSSS.....SSSS.................
    Query:                 .....QQQ..QQQQ...QQQQ
                                |               |
                                |               |
                                queryBinStart   queryBinEnd

    The second example shows how queryStartInSubject can be greater than
    subjectBinStart:

                                queryStartInSubject
                                |               queryEndInSubject
                                |               |
                                |               |
                               subjectBinStart  |subjectBinEnd
                               ||               ||
                               ||               ||
    Subject:  .................SSS..SSSS.....SSSS.................
    Query:                      QQQ..QQQQ...QQQQ
                                |               |
                                |               |
                                queryBinStart   queryBinEnd

    @param bin_: a C{light.histogram.bin}.
    @param queryLen: the length of the query sequence.
    @return: A C{dict} with C{str} keys as follows:

            subjectBinStart: C{int} subject index of the bin start.
            subjectBinEnd: C{int} subject index of the bin end.
            queryBinStart: C{int} query index of the bin start.
            queryBinEnd: C{int} query index of the bin end.
            queryStartInSubject: C{int} subject index of the query start.
            queryEndInSubject: C{int} subject index of the query end.

        The end offsets are given in normal Python style, so for example
        the region of the subject covered by the bin could be extracted from
        the subject via subject[subjectBinStart:subjectBinEnd]
    """

    def debugPrint(locals, msg):
        """
        Print debugging information to stderr showing the local variables from
        a call to normalizeBin and then raise an C{AssertionError}.

        @param locals: A C{dict} of local variables.
        @param msg: A C{str} message to raise C{AssertionError} with.
        @raises: Unconditionally raises an C{AssertionError}.
        """
        print('normalizeBin error:', file=sys.stderr)
        print('  queryLen: %d' % queryLen, file=sys.stderr)
        skipVars = set(('debugPrint', 'bin_',
                        'allQueryFeatures', 'allQueryOffsets',
                        'allSubjectFeatures', 'allSubjectOffsets'))
        for var in sorted(set(locals) - skipVars):
            print('  %s: %s' % (var, locals[var]), file=sys.stderr)

        print('allQueryOffsets:',
              ' '.join(map(str, sorted(allQueryOffsets))), file=sys.stderr)
        print('allSubjectOffsets:',
              ' '.join(map(str, sorted(allSubjectOffsets))), file=sys.stderr)

        # Sort features based only on offset. We can't just use sorted()
        # because the features can be a mix of landmarks and trig points
        # and trig points don't have a symbolDetail attribute, which causes
        # the __lt__ comparison to fail
        key = attrgetter('offset')

        allQueryFeatures = locals['allQueryFeatures']
        if allQueryFeatures:
            print('All query features:', file=sys.stderr)
            for feature in sorted(allQueryFeatures, key=key):
                print('  ', feature, file=sys.stderr)
        else:
            print('All query features is empty.', file=sys.stderr)

        allSubjectFeatures = locals['allSubjectFeatures']
        if allSubjectFeatures:
            print('All subject features:', file=sys.stderr)
            for feature in sorted(allSubjectFeatures, key=key):
                print('  ', feature, file=sys.stderr)
        else:
            print('All subject features is empty.', file=sys.stderr)

        _printBin(bin_)
        raise AssertionError(msg)

    allSubjectFeatures, allSubjectOffsets = histogramBinFeatures(bin_,
                                                                 'subject')
    allQueryFeatures, allQueryOffsets = histogramBinFeatures(bin_, 'query')

    # The following 4 variables all have offset values, which can be used
    # as regular Python start, end values for delimiting regions in
    # strings.
    queryBinStart = min(allQueryOffsets)
    queryBinEnd = max(allQueryOffsets)
    subjectBinStart = min(allSubjectOffsets)
    subjectBinEnd = max(allSubjectOffsets)

    queryStartInSubject = subjectBinStart - queryBinStart
    queryEndInSubject = queryStartInSubject + queryLen

    return {
        'subjectBinStart': subjectBinStart,
        'subjectBinEnd': subjectBinEnd,
        'queryBinStart': queryBinStart,
        'queryBinEnd': queryBinEnd,
        'queryStartInSubject': queryStartInSubject,
        'queryEndInSubject': queryEndInSubject,
    }
