# NOTE: This code is not yet in use. But there are tests for it based on
#       light matter mini-HSP-like objects (the info found in the 'matches'
#       info dicts passed to a Result by Database.find). So it will become
#       useful when we want to display those details.

from light.bin_score import histogramBinFeatures


def printBin(bin_, indent=''):
    for key, value in bin_.items():
        print('%s%s: %s' % (indent, key, value))


def normalizeBin(bin_, queryLen):
    """
    Examine the sense of a bin_ and return information about where the
    query and the bin begin and end.  Return a dict with keys that allow the
    query and the bin to be displayed relative to the subject orientation
    (i.e., with start < stop for both the read and the match). The returned
    query indices are offsets into the subject. I.e., they indicate where on
    the subject the query lies.

    The diagram below shows a possible way in which the query and the subject
    can align and the terminology used. Note that this doesn't display all
    possibilities.

                           queryStartInSubject              queryEndInSubject
                           |                               |
                           |                               |
                           |   subjectBinStart  subjectBinEnd
                           |   |               |           |
                           |   |               |           |
    Subject:  .................MMMMMMMMMMMMMMMM.................
    Query:                 ....MMMMMMMMMMMMMMMM............
                               |               |
                               |               |
                               queryBinStart   queryBinEnd

    NOTE: the returned queryStartInSubject value may be negative.  We
    consider the subject sequence to start at offset 0.  So if the query
    string has sufficient additional residues before the start of the
    alignment match, it may protrude to the left of the subject. Similarly,
    the returned queryEndInSubject can be greater than the subjectEnd.

    @param bin_: a C{light.histogram.bin}.
    @param queryLen: the length of the query sequence.
    """

    def debugPrint(locals, msg=None):
        """
        Print debugging information showing the local variables from
        a call to normalizeHSP and then raise an C{AssertionError}.

        @param locals: A C{dict} of local variables.
        @param msg: A C{str} message to raise C{AssertionError} with.
        """
        print('normalizeBin error:')
        print('  queryLen: %d' % queryLen)
        for var in sorted(locals.keys()):
            if var in ('debugPrint', 'bin'):
                continue
            print('  %s: %s' % (var, locals[var]))
        print('  Original Bin:')
        printBin(bin_, '    ')
        if msg:
            raise AssertionError(msg)
        else:
            raise AssertionError()

    allSubjectFeatures, allSubjectOffsets = histogramBinFeatures(bin_,
                                                                 'subject')
    allQueryFeatures, allQueryOffsets = histogramBinFeatures(bin_, 'query')

    queryBinStart = min(allQueryOffsets)
    queryBinEnd = max(allQueryOffsets)
    subjectBinStart = min(allSubjectOffsets)
    subjectBinEnd = max(allSubjectOffsets)
    queryStartInSubject = subjectBinStart - queryBinStart
    queryEndInSubject = queryStartInSubject + queryLen

    # Sanity checks.
    if queryStartInSubject > subjectBinStart:
        debugPrint(locals(), 'queryStartInSubject > subjectBinStart')
    if queryEndInSubject < subjectBinEnd:
        debugPrint(locals(), 'queryEndInSubject < subjectBinEnd')

    return {
        'subjectBinStart': subjectBinStart,
        'subjectBinEnd': subjectBinEnd,
        'queryBinStart': queryBinStart,
        'queryBinEnd': queryBinEnd,
        'queryStartInSubject': queryStartInSubject,
        'queryEndInSubject': queryEndInSubject,
    }
