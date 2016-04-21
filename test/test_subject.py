import six
from unittest import TestCase
from six import StringIO

from dark.reads import AARead, SSAARead

from light.subject import Subject, SubjectStore
from light.exceptions import SubjectStoreException


class TestSubject(TestCase):
    """
    Tests for the light.database.Subject class.
    """
    def testHashCountIsStored(self):
        """
        A Subject must store the hashcount it is passed.
        """
        read = AARead('id', 'AA')
        self.assertEqual(6, Subject(read, 6).hashCount)

    def testHashCountDefaultsToZero(self):
        """
        A Subject must store a hashcount of zero if the hashcount isn't given.
        """
        read = AARead('id', 'AA')
        self.assertEqual(0, Subject(read).hashCount)

    def testLen(self):
        """
        A Subject must return the correct (sequence) length via len().
        """
        self.assertEqual(2, len(Subject(AARead('id', 'AA'))))

    def testIdenticalSubjectsCompareEqual(self):
        """
        Two identical Subject instances must compare equal.
        """
        read = AARead('id', 'AA')
        self.assertEqual(Subject(read), Subject(read))

    def testSubjectsCompareEqualEvenWithDifferingHashCounts(self):
        """
        Two otherwise identical Subject instances must compare equal even if
        their hash counts differ.
        """
        read = AARead('id', 'AA')
        self.assertEqual(Subject(read, 1), Subject(read, 2))

    def testSubjectsOfIncompatibleReadTypesCompareUnequal(self):
        """
        Two Subject instances that have reads that cannot be compared must
        compare unequal.
        """
        read1 = AARead('id', 'AA')
        read2 = SSAARead('id', 'AA', 'HH')
        self.assertNotEqual(Subject(read1), Subject(read2))


class TestSubjectStore(TestCase):
    """
    Tests for the light.database.SubjectStore class.
    """
    def testEmptyStoreHasZeroLength(self):
        """
        An empty SubjectStore must have zero length.
        """
        self.assertEqual(0, len(SubjectStore()))

    def testLengthIsOneAfterAddSubject(self):
        """
        Adding a subject to a new SubjectStore results in SubjectStore with
        length one.
        """
        ss = SubjectStore()
        subject = Subject(AARead('id', 'AA'))
        ss.add(subject)
        self.assertEqual(1, len(ss))

    def testAddSubjectWithoutIndex(self):
        """
        Adding one subject to a SubjectStore but not giving a subject index
        must return the expected result.
        """
        subject = Subject(AARead('id', 'AA'))
        preExisting, subjectIndex = SubjectStore().add(subject)
        self.assertFalse(preExisting)
        self.assertEqual('0', subjectIndex)

    def testAddSubjectWithIndex(self):
        """
        Adding one subject to a SubjectStore and giving a subject index
        must return the expected result.
        """
        subject = Subject(AARead('id', 'AA'))
        preExisting, subjectIndex = SubjectStore().add(subject, '3')
        self.assertFalse(preExisting)
        self.assertEqual('3', subjectIndex)

    def testAddTwoSubjectsWithIndex(self):
        """
        Adding two subjects to a SubjectStore and giving their indices
        must have the expected result.
        """
        ss = SubjectStore()
        subject1 = Subject(AARead('id1', 'AA'))
        subject2 = Subject(AARead('id2', 'AA'))
        ss.add(subject1, '1')
        ss.add(subject2, '2')
        self.assertEqual(2, len(ss))

    def testAddSubjectTwiceWithSameIndex(self):
        """
        Adding a subject to a SubjectStore and giving a subject index
        and then re-adding it with the same index must return the expected
        result.
        """
        subject = Subject(AARead('id', 'AA'))
        ss = SubjectStore()
        ss.add(subject, '3')
        preExisting, subjectIndex = ss.add(subject, '3')
        self.assertTrue(preExisting)
        self.assertEqual('3', subjectIndex)

    def testLengthIsOneAfterAddSameSubjectTwice(self):
        """
        Adding a subject to a new SubjectStore twice (with an identical index)
        results in SubjectStore with length one.
        """
        ss = SubjectStore()
        subject = Subject(AARead('id', 'AA'))
        ss.add(subject)
        ss.add(subject)
        self.assertEqual(1, len(ss))

    def testAddSubjectTwiceWithDifferentIndex(self):
        """
        Adding a subject to a SubjectStore and giving a subject index
        and then re-adding it with a different index must raise a
        SubjectStoreException.
        """
        subject = Subject(AARead('id', 'AA'))
        ss = SubjectStore()
        ss.add(subject, '3')
        error = ("^Already stored index \(3\) for subject 'id' does not match "
                 "subsequently passed index \(4\) for the subject\.$")
        six.assertRaisesRegex(self, SubjectStoreException, error, ss.add,
                              subject, '4')

    def testGetSubjects(self):
        """
        getSubjects must return the expected result.
        """
        ss = SubjectStore()
        subject1 = Subject(AARead('id1', 'AA'))
        subject2 = Subject(AARead('id2', 'CC'))
        ss.add(subject1, '1')
        ss.add(subject2, '2')
        result = set(ss.getSubjects())
        self.assertEqual(set([subject1, subject2]), result)

    def testSaveRestoreEmpty(self):
        """
        After a save/restore of an empty subject store, the result must be as
        expected.
        """
        ss = SubjectStore()
        fp = StringIO()
        ss.save(fp)
        fp.seek(0)
        ss = SubjectStore.restore(fp)
        self.assertEqual(0, len(ss))

    def testSaveRestoreNonEmpty(self):
        """
        After a save/restore of a non-empty subject store, the subject store
        must be as expected.
        """
        ss = SubjectStore()
        subject = Subject(AARead('id', 'AA'), 6)
        ss.add(subject, '3')
        fp = StringIO()
        ss.save(fp)
        fp.seek(0)
        ss = SubjectStore.restore(fp)
        self.assertEqual(1, len(ss))
        self.assertEqual('3', ss.getIndexBySubject(subject))
        result = ss.getSubjectByIndex('3')
        self.assertEqual(subject, result)
        self.assertEqual(6, result.hashCount)
