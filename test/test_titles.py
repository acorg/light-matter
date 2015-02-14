from unittest import TestCase, skip
from mock import patch

from mocking import mockOpen
from sample_data import (
    DB, PARAMS, COWPOX, MONKEYPOX, MUMMYPOX, SQUIRRELPOX55, SQUIRRELPOX1296,
    RECORD0, RECORD1, RECORD2, RECORD3, RECORD4,
    READ0, READ2, READ3)

from dark.hsp import HSP
from dark.titles import titleCounts, TitleAlignments, TitlesAlignments

from light.alignments import LightReadsAlignments


class TestTitleCounts(TestCase):
    """
    Test the titleCounts function.
    """

    def testEmpty(self):
        """
        If passed an empty readsAlignments, titleCounts must return an
        empty dictionary.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual({}, titleCounts(readsAlignments))

    def testThreeRecords(self):
        """
        If alignments for three reads are passed to titleCounts, it must
        return the correct title counts.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual(
                {
                    SQUIRRELPOX1296.id: 1,
                    MUMMYPOX.id: 1,
                    COWPOX.id: 1,
                    MONKEYPOX.id: 1,
                    SQUIRRELPOX55.id: 1
                },
                titleCounts(readsAlignments))

    def testDuplicatedTitle(self):
        """
        If alignments for reads have a common title, the count on that title
        must be correct.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD2 + RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual(
                {
                    COWPOX.id: 2,
                },
                titleCounts(readsAlignments))


class TestTitlesAlignments(TestCase):
    """
    Test the TitlesAlignments class
    """

    def testEmpty(self):
        """
        An instance of TitlesAlignments must have no titles if passed an
        empty readsAlignments instance.
        """
        mockOpener = mockOpen(read_data=(PARAMS))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertEqual([], titlesAlignments.keys())

    def testExpectedTitles(self):
        """
        An instance of TitlesAlignments must have the expected titles.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertEqual(
                [
                    COWPOX.id,
                    MONKEYPOX.id,
                    MUMMYPOX.id,
                    SQUIRRELPOX1296.id,
                    SQUIRRELPOX55.id,
                ],
                sorted(titlesAlignments.keys()))

    def testExpectedTitleDetails(self):
        """
        An instance of TitleAlignments in a TitlesAlignments instance must
        have the expected attributes.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0))
        with patch('__builtin__.open', mockOpener, create=True):

            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)

            titleAlignments = titlesAlignments[SQUIRRELPOX1296.id]
            self.assertEqual(SQUIRRELPOX1296.id, titleAlignments.subjectTitle)
            self.assertEqual(len(SQUIRRELPOX1296.sequence),
                             titleAlignments.subjectLength)
            self.assertEqual(1, len(titleAlignments))
            self.assertEqual(READ0, titleAlignments[0].read)
            self.assertEqual(HSP(5), titleAlignments[0].hsps[0])

            titleAlignments = titlesAlignments[SQUIRRELPOX55.id]
            self.assertEqual(SQUIRRELPOX55.id, titleAlignments.subjectTitle)
            self.assertEqual(len(SQUIRRELPOX55.sequence),
                             titleAlignments.subjectLength)
            self.assertEqual(1, len(titleAlignments))
            self.assertEqual(READ0, titleAlignments[0].read)
            self.assertEqual(HSP(4), titleAlignments[0].hsps[0])

    def testTitleCollection(self):
        """
        A title that occurs in the alignments of multiple reads must have
        the data from those reads collected properly.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD2 + RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)

            title = COWPOX.id
            titleAlignments = titlesAlignments[title]
            self.assertEqual(2, len(titleAlignments))

            self.assertEqual(title, titleAlignments.subjectTitle)
            self.assertEqual(len(COWPOX.sequence),
                             titleAlignments.subjectLength)

            self.assertEqual(READ2, titleAlignments[0].read)
            self.assertEqual(HSP(10), titleAlignments[0].hsps[0])

            self.assertEqual(READ3, titleAlignments[1].read)
            self.assertEqual(HSP(10), titleAlignments[1].hsps[0])

    def testAddTitleRepeat(self):
        """
        The addTitle function must raise a C{KeyError} if an attempt is made
        to add a pre-existing title to a TitlesAlignments instance.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = SQUIRRELPOX1296.id
            titleAlignments = TitleAlignments(title, 55)
            error = ("Title 'Squirrelpox virus 1296/99' already present in "
                     "TitlesAlignments instance\.")
            self.assertRaisesRegexp(KeyError, error, titlesAlignments.addTitle,
                                    title, titleAlignments)

    def testAddTitle(self):
        """
        The addTitle function must add a title to the TitlesAlignments
        instance.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = 'Squirrelpox virus 23'
            titleAlignments = TitleAlignments(title, 55)
            self.assertTrue(title not in titlesAlignments)
            titlesAlignments.addTitle(title, titleAlignments)
            self.assertTrue(title in titlesAlignments)

    def testHsps(self):
        """
        The hsps function must yield all the hsps for all titles in a
        TitlesAlignments instance.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = list(titlesAlignments.hsps())
            self.assertEqual(
                sorted([HSP(5), HSP(4), HSP(3), HSP(1), HSP(10)]),
                sorted(result))


class TestTitlesAlignmentsFiltering(TestCase):
    """
    Test the TitlesAlignments class filter function.
    """

    def testFilterWithNoArguments(self):
        """
        The filter function must return a TitlesAlignments instance with all
        the titles of the original when called with no arguments.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter()
            self.assertEqual(
                [
                    COWPOX.id,
                    MONKEYPOX.id,
                    MUMMYPOX.id,
                    SQUIRRELPOX1296.id,
                    SQUIRRELPOX55.id,
                ],
                sorted(result.keys()))

    def testMinMatchingReads(self):
        """
        The filter function work correctly when passed a value for
        minMatchingReads.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minMatchingReads=2)
            self.assertEqual(
                [
                    COWPOX.id,
                ],
                result.keys())

    def testMinMedianScore(self):
        """
        The filter function work correctly when passed a value for
        minMedianScore.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minMedianScore=10)
            self.assertEqual(
                [
                    COWPOX.id,
                ],
                result.keys())

    def testWithScoreBetterThan(self):
        """
        The filter function work correctly when passed a value for
        withScoreBetterThan.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(withScoreBetterThan=5)
            self.assertEqual(
                [
                    COWPOX.id,
                ],
                result.keys())

    def testReadSetFilterAllowAnything(self):
        """
        The filter function work correctly when passed a 0.0 value for
        minNewReads, i.e. that considers any read set sufficiently novel.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minNewReads=0.0)
            self.assertEqual(
                [
                    COWPOX.id,
                    MONKEYPOX.id,
                    MUMMYPOX.id,
                    SQUIRRELPOX1296.id,
                    SQUIRRELPOX55.id,
                ],
                sorted(result.keys()))

    def testReadSetFilterStrict(self):
        """
        The filter function work correctly when passed a 1.0 value for
        minNewReads.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minNewReads=1.0)

            # Either MUMMYPOX.id invalidates MONKEYPOX.id or vice-versa. It
            # depends on Python's dict walking order. Check for both,
            # making sure just one of them is true.

            assertionCount = 0
            if MUMMYPOX.id in result:
                self.assertTrue(MONKEYPOX.id in
                                result.readSetFilter.invalidates(MUMMYPOX.id))
                assertionCount += 1
            if MONKEYPOX.id in result:
                self.assertTrue(MUMMYPOX.id in
                                result.readSetFilter.invalidates(MONKEYPOX.id))
                assertionCount += 1

            self.assertEqual(1, assertionCount)

    @skip('Coverage is not yet implemented for light matter title alignments.')
    def testCoverageExcludesAll(self):
        """
        The coverage function must return an titlesAlignments instance with
        no titles if none of its titles has sufficient coverage.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minCoverage=0.1)
            self.assertEqual(0, len(result))

    @skip('Coverage is not yet implemented for light matter title alignments.')
    def testCoverageIncludesAll(self):
        """
        The coverage function must return an titlesAlignments instance with
        all titles if all its titles has sufficient coverage.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minCoverage=0.0)
            self.assertEqual(
                [
                    COWPOX.id,
                    MONKEYPOX.id,
                    MUMMYPOX.id,
                    SQUIRRELPOX1296.id,
                    SQUIRRELPOX55.id,
                ],
                sorted(result.keys()))

    @skip('Coverage is not yet implemented for light matter title alignments.')
    def testCoverageIncludesSome(self):
        """
        The coverage function must return an titlesAlignments instance with
        only the expected titles if only some of its titles have sufficient
        coverage.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            # To understand why the following produces the result it does,
            # you need to look at the HSP coverage in sample_data.py and
            # calculate the coverage by hand.
            result = titlesAlignments.filter(minCoverage=0.0011)
            self.assertEqual(
                [
                    COWPOX.id,
                    MONKEYPOX.id,
                    MUMMYPOX.id,
                ],
                sorted(result.keys()))


class TestTitleSorting(TestCase):
    """
    Tests for the L{dark.titles.TitlesAlignments.sortTitles} function.
    """

    def testUnknown(self):
        """
        Sorting on an unknown attribute must raise C{ValueError}.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertRaises(ValueError, titlesAlignments.sortTitles, 'xxx')

    def testEmpty(self):
        """
        Sorting when there are no titles must return the empty list.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('title')
            self.assertEqual([], result)

    def testMedianScore(self):
        """
        Sorting on median score must work, including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3 + RECORD4))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('medianScore')
            self.assertEqual([
                COWPOX.id,           # 24
                SQUIRRELPOX1296.id,  # 2
                SQUIRRELPOX55.id,    # 2
                MONKEYPOX.id,        # 1
                MUMMYPOX.id,         # 1
            ], result)

    def testMaxScore(self):
        """
        Sorting on max score must work, including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('maxScore')
            self.assertEqual([
                COWPOX.id,           # 24
                SQUIRRELPOX1296.id,  # 2
                SQUIRRELPOX55.id,    # 2
                MONKEYPOX.id,        # 1
                MUMMYPOX.id,         # 1
            ], result)

    def testReadCount(self):
        """
        Sorting on read count must work, including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('readCount')
            self.assertEqual([
                COWPOX.id,           # 3
                MONKEYPOX.id,        # 1
                MUMMYPOX.id,         # 1
                SQUIRRELPOX1296.id,  # 1
                SQUIRRELPOX55.id,    # 1
            ], result)

    def testLength(self):
        """
        Sorting on sequence length must work.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('length')
            self.assertEqual([
                COWPOX.id,           # 21
                SQUIRRELPOX1296.id,  # 15
                MONKEYPOX.id,        # 13
                SQUIRRELPOX55.id,    # 11
                MUMMYPOX.id,         # 8
            ], result)

    def testTitle(self):
        """
        Sorting on title must work.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch('__builtin__.open', mockOpener, create=True):
            readsAlignments = LightReadsAlignments('file.json', DB)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('title')
            self.assertEqual([
                COWPOX.id,
                MONKEYPOX.id,
                MUMMYPOX.id,
                SQUIRRELPOX1296.id,
                SQUIRRELPOX55.id,
            ], result)
