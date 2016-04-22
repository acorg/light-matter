import sys
try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from dark.reads import readClassNameToClass

from light.exceptions import SubjectStoreException


class Subject(object):
    """
    Hold information about a database subject.

    @param read: A C{Read} instance, or an instance of a C{Read} subclass.
    @param hashCount: An C{int} count of the number of hashes found in the
        subject when it was added (via Database.addSubject).
    """
    def __init__(self, read, hashCount=0):
        self.read = read
        self.hashCount = hashCount

    def __eq__(self, other):
        """
        Are two subjects equal?

        @return: A C{bool} to indicate equality.
        """
        try:
            return self.read == other.read
        except AttributeError:
            # This occurs when two subjects have read types that cannot be
            # compared (e.g., an SSAARead has a structure but no quality,
            # and an AARead has quality but no structure). In normal
            # operation this will not happen (all reads are likely to be
            # AARead or AAReadWithX instances), but there's no reason reads
            # of non-comparable types can't be added to a database if
            # someone wants to do that during development.
            return False

    def __hash__(self):
        """
        Calculate a hash key for a subject.

        @return: The C{int} hash key for the subject.
        """
        return self.read.__hash__()

    def __len__(self):
        """
        Get a subject's length.

        @return: The C{int} length of the subject's read sequence.
        """
        return len(self.read)


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
                        storedIndex, subject.read.id, subjectIndex))
            return (True, storedIndex)

    def getSubjectByIndex(self, subjectIndex):
        """
        Return information about a subject, given its index.

        @param subjectIndex: A C{str} subject index.
        @raises KeyError: if the subject index is unknown (i.e., is not present
            in C{self._indexToSubject}).
        @return: A C{Subject} instance.
        """
        return self._indexToSubject[subjectIndex]

    def getIndexBySubject(self, subject):
        """
        Return a subject's index.

        @param subject: A C{Subject} instance whose index is wanted.
        @raises KeyError: if the passed C{Subject} instance is unknown (i.e.,
            is not present in C{self._subjectToIndex}).
        @return: A C{str} subject index.
        """
        try:
            return self._subjectToIndex[subject]
        except KeyError:
            # Be user friendly and raise a key error containing the subject
            # id, instead of the md5 hash value of the subject.
            raise KeyError(subject.read.id)

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
        state = []
        for subject, subjectIndex in self._subjectToIndex.items():
            state.append({
                'readClass': subject.read.__class__.__name__,
                'hashCount': subject.hashCount,
                'readDict': subject.read.toDict(),
                'subjectIndex': subjectIndex,
            })

        dump(state, fp)
        fp.write('\n')

    @classmethod
    def restore(cls, fp=sys.stdin):
        """
        Load a subject store from a file.

        @param fp: A file pointer.
        @raises ValueError: If valid JSON cannot be loaded from C{fp}.
        @raises KeyError: If the subject store contains an unknown read class.
        @raises SubjectStoreException: If a restored subject was already
            present in the store or if a subject cannot be stored with its
            original subject index. Both of which should be impossible.
        @return: A new instance of L{SubjectStore}.
        """
        new = cls()
        add = new.add

        state = loads(fp.readline()[:-1])

        for subjectInfo in state:
            readClass = readClassNameToClass[subjectInfo['readClass']]
            read = readClass.fromDict(subjectInfo['readDict'])
            subject = Subject(read, subjectInfo['hashCount'])
            subjectIndex = subjectInfo['subjectIndex']
            preExisting, actualSubjectIndex = add(subject, subjectIndex)

            # Sanity checks.
            if preExisting:
                raise SubjectStoreException(
                    'Subject id %r already present in subject store during '
                    'restore!' % read.id)

            if actualSubjectIndex != subjectIndex:
                raise SubjectStoreException(
                    'Attempted to add subject id %r to subject store with '
                    'subject index %r during restore, but the subject was '
                    'actually stored with a different index, %r.' % (
                        read.id, subjectIndex, actualSubjectIndex))

        return new
