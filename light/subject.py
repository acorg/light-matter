import sys
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from dark.reads import AARead

from light.exceptions import SubjectStoreException


class Subject(AARead):
    """
    Hold information about a database subject.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information (might be
        nucleotides or proteins).
    @param hashCount: An C{int} count of the number of hashes found in the
        subject when it was added (via Database.addSubject).
    @param quality: An optional C{str} of phred quality scores. If not C{None},
        it must be the same length as C{sequence}.
    """
    def __init__(self, id, sequence, hashCount, quality=None):
        AARead.__init__(self, id, sequence, quality)
        self.hashCount = hashCount


class SubjectStore(object):
    """
    Maintain information about the set of subjects in a database.
    """
    def __init__(self):
        self._indexToSubject = {}  # maps str indices to subjects.
        self._subjectToIndex = {}  # maps subjects to str indices.

    def __len__(self):
        return len(self._subjectToIndex)

    def add(self, subject, subjectIndex=None):
        """
        Add a subject to our collection.

        @param subject: A C{Subject} instance.
        @param subjectIndex: A C{str} index to use for this subject. If
            C{None} an index should be chosen.
        @return: A 2-tuple (pre-existing, index) where 'pre-existing' is a
            C{bool} to indicate whether the subject was already present in the
            store, and 'index' is the C{str} subject index of the added
            subject.
        """
        try:
            storedIndex = self._subjectToIndex[subject]
        except KeyError:
            if subjectIndex is None:
                subjectIndex = str(len(self))
            self._subjectToIndex[subject] = subjectIndex
            self._indexToSubject[subjectIndex] = subject
            return (False, subjectIndex)
        else:
            if subjectIndex is not None and storedIndex != subjectIndex:
                raise SubjectStoreException(
                    'Already stored index (%s) for subject %r does not '
                    'match subsequently passed index (%s) for the subject.' % (
                        storedIndex, subject.id, subjectIndex))
            return (True, storedIndex)

    def getSubjectByIndex(self, subjectIndex):
        """
        Return information about a subject, given its index.

        @param subjectIndex: A C{str} subject index.
        @return: A C{Subject} instance.
        @raises KeyError: if the subject index is unknown.
        """
        return self._indexToSubject[subjectIndex]

    def getIndexBySubject(self, subject):
        """
        Return a subject's index.

        @param subject: A C{Subject} instance whose index is wanted.
        @return: A C{str} subject index.
        @raises KeyError: if the passed C{Subject} is unknown.
        """
        try:
            return self._subjectToIndex[subject]
        except KeyError:
            # Be user friendly and raise a key error containing the subject
            # id, instead of the md5 sum of the id and sequence.
            raise KeyError(subject.id)

    def getSubjects(self):
        """
        Return information about all subjects.

        @return: a generator that yields C{Subject} instances.
        """
        return self._subjectToIndex.keys()

    def save(self, fp=sys.stdout):
        """
        Save the subject store to a file.

        @param fp: A file pointer.
        """
        state = [[subjectIndex, s.id, s.sequence, s.hashCount, s.quality]
                 for s, subjectIndex in self._subjectToIndex.items()]
        dump(state, fp)
        fp.write('\n')

    @classmethod
    def restore(cls, fp=sys.stdin):
        """
        Load a subject store from a file.

        @param fp: A file pointer.
        @return: An instance of L{SubjectStore}.
        @raises ValueError: If valid JSON cannot be loaded from C{fp}.
        """
        new = cls()
        state = loads(fp.readline()[:-1])

        add = new.add
        for subjectIndex, sequenceId, sequence, hashCount, quality in state:
            add(Subject(sequenceId, sequence, hashCount, quality),
                subjectIndex)

        return new
