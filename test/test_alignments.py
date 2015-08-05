import builtins
from json import dumps, loads
import bz2
from io import StringIO
from unittest import TestCase, skip
from unittest.mock import patch

from .mocking import mockOpen
from .sample_data import (
    DB, PARAMS, RECORD0, RECORD1, RECORD2, RECORD3,
    CATPOX, COWPOX, MUMMYPOX, MONKEYPOX, SQUIRRELPOX,
    READ1, READ2,
    READ0_SQUIRRELPOX_SCORE, READ0_CATPOX_SCORE,
    READ1_MONKEYPOX_SCORE, READ1_MONKEYPOX_HSP2_SCORE, READ1_MUMMYPOX_SCORE,
    READ2_COWPOX_SCORE,
    READ3_COWPOX_SCORE)

from dark.hsp import HSP
from dark.reads import AARead

from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.beta_strand import BetaStrand
from light.trig.amino_acids import AminoAcids
from light.database import Database
from light.parameters import Parameters
from light.alignments import LightReadsAlignments, jsonDictToAlignments


class BZ2(object):
    """
    A BZ2File mock.
    """
    def __init__(self, data):
        self._data = data
        self._index = 0

    def close(self):
        pass

    def readline(self):
        self._index += 1
        return self._data[self._index - 1]

    def __iter__(self):
        index = self._index
        self._index = len(self._data)
        return iter(self._data[index:])


class TestLightReadsAlignments(TestCase):
    """
    Test the LightReadsAlignments class.
    """

    def testEmptyJSONInput(self):
        """
        When a JSON input file is empty, a C{ValueError} must be raised
        on trying to read it.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Could not convert first line of 'file.json' to "
                     "JSON \(Expected object or value\)\.$")
            self.assertRaisesRegex(
                ValueError, error, LightReadsAlignments, 'file.json', DB)

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the light matter results from it must raise a C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not JSON\n')
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Could not convert first line of 'file\.json' to JSON "
                     "\(Unexpected character found when decoding 'null'\)\.$")
            self.assertRaisesRegex(
                ValueError, error, LightReadsAlignments, 'file.json', DB)

    def testScoreTitle(self):
        """
        The score title must be as expected.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual('Score', readsAlignments.params.scoreTitle)

    def testSubjectIsNucleotides(self):
        """
        The subject is nucleotide parameter must be false.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual('light', readsAlignments.params.application)
            self.assertFalse(readsAlignments.params.subjectIsNucleotides)

    def testApplicationParams(self):
        """
        Light matter parameters must be extracted from the input JSON file and
        stored correctly.
        """
        params = Parameters([AlphaHelix, BetaStrand], [AminoAcids])
        savedParams = params.save(StringIO()).getvalue()

        mockOpener = mockOpen(read_data=savedParams)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            diffs = readsAlignments.params.applicationParams.compare(params)
            self.assertIs(diffs, None)

    def testJSONParamsButNoHits(self):
        """
        When light matter parameters are present in the input but there are no
        records, the __iter__ method of a L{LightReadsAlignments} instance must
        not yield anything.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual([], list(readsAlignments))

    def testOneCompressedJSONInput(self):
        """
        If a compressed (bz2) JSON file contains a parameters section and one
        record, it must be read correctly.
        """
        result = BZ2([PARAMS, RECORD0])

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.return_value = result
            readsAlignments = LightReadsAlignments('file.json.bz2', DB)
            self.assertEqual(1, len(list(readsAlignments)))

    def testTwoCompressedJSONInputs(self):
        """
        If two compressed (bz2) JSON files are passed to
        L{LightReadsAlignments} each with a parameters section and one
        record, both records must be read correctly and the result should
        have 2 records.
        """

        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([PARAMS, RECORD0])
                else:
                    return BZ2([PARAMS, RECORD1])

        sideEffect = SideEffect()
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            readsAlignments = LightReadsAlignments(
                ['file1.json.bz2', 'file2.json.bz2'], DB)
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual('read1', result[1].read.id)

    def testThreeCompressedJSONInputs(self):
        """
        If three compressed (bz2) JSON files are passed to
        L{LightReadsAlignments} with names that have a numeric prefix and
        each with a parameters section and one record, all records must be
        read correctly and the result should have 3 records in the correct
        order.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename):
                if self.count == 0:
                    self.test.assertEqual('1.json.bz2', filename)
                    self.count += 1
                    return BZ2([PARAMS, RECORD0])
                elif self.count == 1:
                    self.test.assertEqual('2.json.bz2', filename)
                    self.count += 1
                    return BZ2([PARAMS, RECORD1])
                else:
                    self.test.assertEqual('3.json.bz2', filename)
                    return BZ2([PARAMS, RECORD2])

        sideEffect = SideEffect(self)
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect

            # Note the files are given out of order. Their names will be
            # sorted before they are opened. The sorting of the names is
            # verified in the SideEffect class, above.
            readsAlignments = LightReadsAlignments(
                ['3.json.bz2', '1.json.bz2', '2.json.bz2'], DB)
            result = list(readsAlignments)
            self.assertEqual(3, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual('read1', result[1].read.id)
            self.assertEqual('read2', result[2].read.id)

    def testIncompatibleParameters(self):
        """
        If two compressed (bz2) JSON files with incompatible parameters
        are given to L{LightReadsAlignments}, a C{ValueError} must be
        raised when the files are read.
        """

        class SideEffect(object):
            def __init__(self):
                self._first = True

            def sideEffect(self, _ignoredFilename):
                if self._first:
                    self._first = False
                    return BZ2([PARAMS, RECORD0])
                else:
                    params = loads(PARAMS)
                    params['limitPerLandmark'] = 100
                    return BZ2([dumps(params) + '\n', RECORD1])

        sideEffect = SideEffect()
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            error = ("^Incompatible light matter parameters found\. The "
                     "parameters in file2\.json\.bz2 differ from those "
                     "originally found in file1\.json\.bz2\. Summary of "
                     "differences:\n\tParam 'limitPerLandmark' values 10 "
                     "and 100 differ\.$")
            readsAlignments = LightReadsAlignments(
                ['file1.json.bz2', 'file2.json.bz2'], DB)
            self.assertRaisesRegex(ValueError, error, list, readsAlignments)

    def testGetSubjectSequence(self):
        """
        The getSubjectSequence function must return a correct C{SeqIO.read}
        instance.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            subject = readsAlignments.getSubjectSequence(COWPOX.id)
            self.assertEqual(COWPOX.sequence, subject.sequence)
            self.assertEqual(COWPOX.id, subject.id)

    def testHsps(self):
        """
        The hsps function must yield the HSPs.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2 +
                                         RECORD3))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            self.assertEqual(
                sorted([HSP(READ0_SQUIRRELPOX_SCORE),
                        HSP(READ0_CATPOX_SCORE),
                        HSP(READ1_MONKEYPOX_SCORE),
                        HSP(READ1_MONKEYPOX_HSP2_SCORE),
                        HSP(READ1_MUMMYPOX_SCORE),
                        HSP(READ2_COWPOX_SCORE),
                        HSP(READ3_COWPOX_SCORE)]),
                sorted(readsAlignments.hsps()))


class TestLightReadsAlignmentsFiltering(TestCase):
    """
    Test the LightReadsAlignments class filter function.
    """

    def testNoResultNoFilteringArgs(self):
        """
        If the L{LightReadsAlignments} filter function is called with no
        arguments, and there are no hits, it should produce a generator
        that yields no result.
        """
        mockOpener = mockOpen(read_data=PARAMS)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter())
            self.assertEqual(0, len(result))

    def testOneHitNoFilteringArgs(self):
        """
        If the L{LightReadsAlignments} filter function is called with no
        arguments, and there is one hit, it should produce a generator that
        yields that hit.
        """
        mockOpener = mockOpen(read_data=PARAMS + RECORD0)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter())
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)

    def testLimitZero(self):
        """
        If L{LightReadsAlignments} is limited to zero results, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=PARAMS + RECORD0)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(limit=0))
            self.assertEqual(0, len(result))

    def testLimitOne(self):
        """
        If L{LightReadsAlignments} is limited to one hit, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=PARAMS + RECORD0 + RECORD1)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(limit=1))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)

    def testOneAlignmentPerRead(self):
        """
        If L{LightReadsAlignments} is asked to deliver only the best alignment
        for each read, that must be respected.
        """
        mockOpener = mockOpen(read_data=PARAMS + RECORD0)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(oneAlignmentPerRead=True))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual(SQUIRRELPOX.id, result[0][0].subjectTitle)

    def testScoreCutoffRemovesEntireAlignment(self):
        """
        If the L{LightReadsAlignments} filter function is supposed to filter on
        a scoreCutoff and the cut-off value results in an alignment with no
        HSPs, then the alignment must be removed entirely.
        """
        mockOpener = mockOpen(read_data=PARAMS + RECORD0)
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(
                scoreCutoff=READ0_SQUIRRELPOX_SCORE - 0.01))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual(SQUIRRELPOX.id, result[0][0].subjectTitle)

    def testTitleByRegexCaseInvariant(self):
        """
        Filtering with a title regex must work independent of case.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(titleRegex='sqUIRRel'))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual(SQUIRRELPOX.id, result[0][0].subjectTitle)

    def testTitleByRegexAllAlignments(self):
        """
        Filtering with a title regex must work in the case that all alignments
        for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(titleRegex='squirrel'))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual(SQUIRRELPOX.id, result[0][0].subjectTitle)

    def testTitleByRegexOneAlignments(self):
        """
        Filtering with a title regex must work in the case that only some
        alignments for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(titleRegex='Mummy'))
            self.assertEqual(1, len(result))
            self.assertEqual('read1', result[0].read.id)
            self.assertEqual(MUMMYPOX.id, result[0][0].subjectTitle)

    def testTitleByNegativeRegexOneAlignment(self):
        """
        Filtering with a negative title regex must work in the case that only
        some alignments for a hit are ruled out (in which case only those
        alignments must be removed but the hit is still valid).
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(negativeTitleRegex='Mummy'))
            self.assertEqual(3, len(result))
            self.assertEqual('read1', result[1].read.id)
            self.assertEqual(1, len(result[1]))
            self.assertEqual(MONKEYPOX.id, result[1][0].subjectTitle)

    def testTitleByNegativeRegexMatchesAll(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and return an empty result.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(negativeTitleRegex='pox'))
            self.assertEqual(0, len(result))

    def testTitleByNegativeRegexMatchingAllWithWhitelist(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and result in no hits, except for any
        whitelisted titles.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            title = SQUIRRELPOX.id
            result = list(readsAlignments.filter(negativeTitleRegex='pox',
                                                 whitelist=[title]))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(title, result[0][0].subjectTitle)

    def testTitleByRegexMatchingAllWithBlacklist(self):
        """
        Filtering with a title regex that matches all alignments
        must keep everything, except for any blacklisted titles.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            blacklist = [SQUIRRELPOX.id, CATPOX.id]
            result = list(readsAlignments.filter(titleRegex='pox',
                                                 blacklist=blacklist))
            self.assertEqual(2, len(result))
            self.assertEqual('read1', result[0].read.id)
            self.assertEqual('read2', result[1].read.id)

    def testMinTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(
                minSequenceLen=len(MUMMYPOX.sequence)))
            self.assertEqual(1, len(result))
            self.assertEqual(READ1.id, result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(MUMMYPOX.id, result[0][0].subjectTitle)

    def testMinTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length and if nothing sufficiently long matches, an empty list of
        alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(minSequenceLen=1000000))
            self.assertEqual(0, len(result))

    def testMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(
                maxSequenceLen=len(COWPOX.sequence) + 1))
            self.assertEqual(1, len(result))
            self.assertEqual(READ2.id, result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(COWPOX.id, result[0][0].subjectTitle)

    def testMaxTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length and if no sufficiently short sequences match, an empty
        list of alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(
                maxSequenceLen=len(COWPOX.sequence) - 1))
            self.assertEqual(0, len(result))

    def testMinAndMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments simultaneously on minimum and
        maximum hit sequence length.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(
                minSequenceLen=len(COWPOX.sequence),
                maxSequenceLen=len(COWPOX.sequence)))
            self.assertEqual(1, len(result))
            self.assertEqual(READ2.id, result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(COWPOX.id, result[0][0].subjectTitle)

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testMinStart(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(minStart=15300))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testMinStartNoHits(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence, and if no hsps match then an empty result set
        must be returned.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(minStart=100000))
            self.assertEqual(0, len(result))

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testMaxStop(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(maxStop=1500))
            self.assertEqual(1, len(result))
            self.assertEqual('read2', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(COWPOX.id, result[0][0].subjectTitle)

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testMaxStopNoHits(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence, and if no hsps match then an empty result set must
        be returned.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(maxStop=100))
            self.assertEqual(0, len(result))

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testMinStartAndMaxstop(self):
        """
        It must be possible to filter alignments based simultaneously on
        mininum and maximum offset in the hit sequence.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(minStart=9000, maxStop=12000))
            self.assertEqual(1, len(result))
            self.assertEqual('read1', result[0].read.id)
            self.assertEqual(2, len(result[0]))

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testRepeatedFilter_MinStartThenMinStart(self):
        """
        It must be possible to filter alignments multiple times using the same
        filter parameters.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            readsAlignments.filter(minStart=9000)
            readsAlignments.filter(minStart=9000)
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('read0', result[0].read.id)
            self.assertEqual('read1', result[1].read.id)

    @skip('Subject filtering on minStart/maxStop not yet implemented.')
    def testRepeatedFilter_MinStartThenMaxstop(self):
        """
        It must be possible to filter alignments multiple times using different
        filter parameters.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            readsAlignments.filter(minStart=9000)
            readsAlignments.filter(maxStop=12000)
            result = list(readsAlignments)
            self.assertEqual(1, len(result))
            self.assertEqual('read1', result[0].read.id)
            self.assertEqual(2, len(result[0]))

    def testClearFilter(self):
        """
        It must be possible to clear any filtering that has been applied.
        """
        result = lambda a: BZ2([PARAMS, RECORD0, RECORD1, RECORD2])

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = result
            readsAlignments = LightReadsAlignments('file.json.bz2', DB)
            self.assertEqual(3, len(list(readsAlignments)))
            readsAlignments.filter(minSequenceLen=14)
            readsAlignments.filter(maxSequenceLen=16)
            readsAlignments.filter(scoreCutoff=0.05)
            result = list(readsAlignments)
            self.assertEqual(1, len(result))
            readsAlignments.clearFilter()
            self.assertEqual(3, len(list(readsAlignments)))

    def testReadIdNoMatches(self):
        """
        When filtering on alignments based on a regex for
        read ids that matches no ids, an empty generator must be returned.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(readIdRegex='blah'))
            self.assertEqual(0, len(result))

    def testReadId(self):
        """
        It must be possible to filter alignments based on a regex for
        read ids.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(readIdRegex='read[12]'))
            self.assertEqual(2, len(result))
            self.assertEqual('read1', result[0].read.id)
            self.assertEqual('read2', result[1].read.id)

    def testReadIdAnchored(self):
        """
        It must be possible to filter alignments based on a regex for
        read ids that is anchored at start and end.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(readIdRegex='^read0$'))
            self.assertEqual(1, len(result))
            self.assertEqual('read0', result[0].read.id)

    def testReadIdCaseSensitive(self):
        """
        Filtering alignments based on a regex for read ids must be case
        sensitive.
        """
        mockOpener = mockOpen(read_data=(PARAMS + RECORD0 + RECORD1 + RECORD2))
        with patch.object(builtins, 'open', mockOpener):
            readsAlignments = LightReadsAlignments('file.json', DB)
            result = list(readsAlignments.filter(readIdRegex='^READ0$'))
            self.assertEqual(0, len(result))


class TestJSONDictToAlignments(TestCase):
    """
    Test the jsonDictToAlignments function.
    """

    def testNoAlignments(self):
        """
        The jsonDictToAlignments function must return an empty list if passed
        a light matter result dict (as would be read from a JSON file) with
        no alignments.
        """
        params = Parameters([], [])
        database = Database(params)
        lightDict = {
            'alignments': [],
        }
        self.assertEqual([], jsonDictToAlignments(lightDict, database))

    def testTwoAlignments(self):
        """
        The jsonDictToAlignments function must produce the correct result when
        passed a light matter result dict (as would be read from a JSON file)
        that has two alignments in it.
        """
        params = Parameters([], [])
        database = Database(params)
        subject0 = AARead('subject0', 'AAA')
        subject1 = AARead('subject1', 'VVV')
        database.addSubject(subject0)
        database.addSubject(subject1)
        lightDict = {
            'alignments': [
                {
                    'subjectIndex': '0',
                    'matchScore': 1.0,
                    'hsps': [
                        {
                            'hspInfo': [
                                {
                                    'landmark': 'AlphaHelix',
                                    'subjectOffset': 0,
                                    'trigPoint': 'Peaks',
                                    'queryOffset': 0,
                                },
                            ],
                            'score': 1.0,
                        },
                    ],
                },
                {
                    'subjectIndex': '1',
                    'matchScore': 1.0,
                    'hsps': [
                        {
                            'hspInfo': [
                                {
                                    'landmark': 'AlphaHelix',
                                    'subjectOffset': 27,
                                    'trigPoint': 'Peaks',
                                    'queryOffset': 27,
                                },
                            ],
                            'score': 1.0,
                        },
                    ],
                },
            ],
            'queryId': 'id',
            'querySequence': 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF',
        }
        alignments = jsonDictToAlignments(lightDict, database)
        self.assertEqual(2, len(alignments))

        self.assertEqual(HSP(1.0), alignments[0].hsps[0])
        self.assertEqual('subject0', alignments[0].subjectTitle)
        self.assertEqual(lightDict['alignments'][0]['hsps'][0]['hspInfo'],
                         alignments[0].hsps[0].hspInfo)

        self.assertEqual(HSP(1.0), alignments[1].hsps[0])
        self.assertEqual('subject1', alignments[1].subjectTitle)
        self.assertEqual(lightDict['alignments'][1]['hsps'][0]['hspInfo'],
                         alignments[1].hsps[0].hspInfo)
