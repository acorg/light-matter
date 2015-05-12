from dark.reads import AARead


class Subject(AARead):
    """
    Hold information about a database subject.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information (might be
        nucleotides or proteins).
    @param hashCount: An C{int} count of the number of hashes found in the
        subject when it was added (via addSubject).
    @param quality: An optional C{str} of phred quality scores. If not C{None},
        it must be the same length as C{sequence}.
    """
    def __init__(self, id, sequence, hashCount, quality=None):
        AARead.__init__(self, id, sequence, quality)
        self.hashCount = hashCount
