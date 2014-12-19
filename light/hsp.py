def printHSP(hsp, indent=''):
    for key, value in hsp.iteritems():
        print '%s%s: %s' % (indent, key, value)


def normalizeHSP(hsp, readLen):
    """
    Examine the sense of an HSP and return information about where the
    read and the alignment (match) begin and end.  Return a dict with keys
    that allow the read and the alignment to be displayed relative to the
    hit orientation (i.e., with start < stop for both the read and the
    match). The returned read indices are offsets into the hit. I.e.,
    they indicate where on the hit the read lies.

    NOTE: the returned readStartInSubject value may be negative.  We
    consider the hit sequence to start at offset 0.  So if the read string
    has sufficient additional residues before the start of the alignment
    match, it may protrude to the left of the hit. Similarly, the returned
    readEndInSubject can be greater than the subjectEnd.

    @param hsp: an HSP in the form of a C{dict}, built from a light matter
        JSON result.
    @param readLen: the length of the read sequence.
    """

    def debugPrint(locals, msg=None):
        """
        Print debugging information showing the local variables from
        a call to normalizeHSP and then raise an C{AssertionError}.

        @param locals: A C{dict} of local variables.
        @param msg: A C{str} message to raise C{AssertionError} with.
        """
        print 'normalizeHSP error:'
        print '  readLen: %d' % readLen
        for var in sorted(locals.keys()):
            if var in ('debugPrint', 'hsp'):
                continue
            print '  %s: %s' % (var, locals[var])
        print '  Original HSP:'
        printHSP(hsp, '    ')
        if msg:
            raise AssertionError(msg)
        else:
            raise AssertionError()

    readStart = hsp['readOffset']
    readEnd = readStart + hsp['landmarkLength']
    subjectStart = hsp['subjectOffset']
    subjectEnd = subjectStart + hsp['landmarkLength']
    readStartInSubject = subjectStart - readStart
    readEndInSubject = readStartInSubject + readLen

    # Sanity checks.
    if readStartInSubject > subjectStart:
        debugPrint(locals(), 'readStartInSubject > subjectStart')
    if readEndInSubject < subjectEnd:
        debugPrint(locals(), 'readEndInSubject < subjectEnd')

    return {
        'readStart': readStart,
        'readEnd': readEnd,
        'readStartInSubject': readStartInSubject,
        'readEndInSubject': readEndInSubject,
        'subjectStart': subjectStart,
        'subjectEnd': subjectEnd,
    }
