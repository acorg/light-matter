from __future__ import division

import six
from unittest import TestCase

from .binScoreTemplate import Template, Query, Subject

from light.features import Landmark, TrigPoint
from light.landmarks import AlphaHelix, EukaryoticLinearMotif
from light.parameters import DatabaseParameters, FindParameters
from light.trig import AminoAcids


# If you edit the following sample template in any way, tests below will
# probably break.
SAMPLE_TEMPLATE = """

    Subject                            FRRRFRRRF-W------|FRRRFRRRFW---RIKR-W-W--W---|--
    AlphaHelix                         FRRRFRRRF        |FRRRFRRRF                  |
    AminoAcids                                   W      |         W        W W      |
    EukaryoticLinearMotif, AminoAcids                   |             RIKR      W   |
    AlphaHelix, EukaryoticLinearMotif                   |FRRRFRRRF    RIKR          |


    Query                                           S-TM---|NGRSFRRRFRRRFTW-T-S-|FRRRFRRRFS-
    Peaks                                           S      |   S              S |         S
    IndividualPeaks                                   T    |             T  T   |
    IndividualTroughs                                  M   |                    |
    AlphaHelix                                             |                    |FRRRFRRRF
    AlphaHelix, AminoAcids                                 |    FRRRFRRRF W     |
    EukaryoticLinearMotif, AminoAcids                      |NGR           W     |

"""


class TestTemplateToList(TestCase):
    """
    Tests for the static templateToList method of the Template class.
    """

    def testEmptyString(self):
        """
        templateToList must return an empty list when passed an empty string.
        """
        self.assertEqual([], Template.templateToList(''))

    def testOneString(self):
        """
        templateToList must return a list of one string when passed a string
        containing no newlines.
        """
        self.assertEqual(['xxx'], Template.templateToList('xxx'))

    def testTwoStrings(self):
        """
        templateToList must return a list of two strings when passed a string
        containing one newline.
        """
        self.assertEqual(['xxx', 'yyy'], Template.templateToList('xxx\nyyy'))

    def testTwoStringsOneLeadingNewline(self):
        """
        templateToList must return a list of two strings when passed a string
        containing two newlines (one at the very start).
        """
        self.assertEqual(['xxx', 'yyy'], Template.templateToList('\nxxx\nyyy'))

    def testTwoStringsTwoLeadingNewlines(self):
        """
        templateToList must return a list of two strings when passed a string
        containing three newlines (two at the very start).
        """
        self.assertEqual(['xx', 'yy'], Template.templateToList('\n\nxx\nyy'))

    def testTwoStringsOneTrailingNewline(self):
        """
        templateToList must return a list of two strings when passed a string
        containing two newlines (one at the very end).
        """
        self.assertEqual(['xxx', 'yyy'], Template.templateToList('xxx\nyyy\n'))

    def testTwoStringsTwoTrailingNewlines(self):
        """
        templateToList must return a list of two strings when passed a string
        containing three newlines (two at the very end).
        """
        self.assertEqual(['xx', 'yy'], Template.templateToList('xx\nyy\n\n'))

    def testStripping(self):
        """
        templateToList must strip whitespace on the right and preserve it on
        the left.
        """
        self.assertEqual(['  x', ' yy'],
                         Template.templateToList('  x  \n yy\t  '))


class TestQueryOrSubjectMixin(object):
    """
    Common tests for the Query and Subject classes.
    """

    def testEmptyTemplate(self):
        """
        The class must raise ValueError when passed an empty template.
        """
        error = "^Could not find '%s' line in template$" % self.CLASS.TYPE
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, [])

    def testTemplateWithNoStart(self):
        """
        The class must raise ValueError when passed a template that does not
        have a query or subject line that indicates the start of the section
        we're interested in.
        """
        error = "^Could not find '%s' line in template$" % self.CLASS.TYPE
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, [])


class TestQuery(TestQueryOrSubjectMixin, TestCase):
    """
    Tests for the Query class.
    """
    CLASS = Query

    def testAttributes(self):
        """
        The class must have the expected attributes.
        """
        self.assertEqual('query', self.CLASS.TYPE)
        self.assertEqual('subject', self.CLASS.OTHER_TYPE)

    def testIndentLength(self):
        """
        The correct indent length must be found for the query.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(52, match.query.indentLength)

    def testRelevantLines(self):
        """
        The relevant lines from the template must be found for the query.
        """
        match = Template(SAMPLE_TEMPLATE)
        expected = [
            '    Query                                           S-TM---|NGRSFRRRFRRRFTW-T-S-|FRRRFRRRFS-',
            '    Peaks                                           S      |   S              S |         S',
            '    IndividualPeaks                                   T    |             T  T   |',
            '    IndividualTroughs                                  M   |                    |',
            '    AlphaHelix                                             |                    |FRRRFRRRF',
            '    AlphaHelix, AminoAcids                                 |    FRRRFRRRF W     |',
            '    EukaryoticLinearMotif, AminoAcids                      |NGR           W     |',
        ]
        self.assertEqual(expected, match.query.relevantLines)

    def testFindPipesMissing(self):
        """
        If a query line doesn't have 2 pipe characters, a ValueError must be
        raised.
        """
        lines = [
            'Subject FRRRF|RRR|FRRRFW',
            'Query STM|STWTS',
        ]
        error = ("^Template line 'Query STM\|STWTS' in query does not "
                 "contain two pipe characters$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testFindPipesNotAligned(self):
        """
        If the pipe characters in a query line don't all align, a ValueError
        must be raised.
        """
        lines = [
            'Subject STM|ST|WTS',
            'Query    FRR|RF|RRR',
            'AminoAcids   W| |',
        ]
        error = ("^Template line 'AminoAcids   W| |' in query has pipe "
                 "characters at offsets 14 and 16, but these do not match the "
                 "offsets \(14 and 17\) of the pipe characters found on the "
                 "first query template line$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testFindPipes(self):
        """
        The offsets of the pipe characters must be calculated correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(59, match.query.pipe1)
        self.assertEqual(80, match.query.pipe2)

    def testNoPipes(self):
        """
        The relevant lines from the template must be found for the query and
        their pipe characters must be removed.
        """
        match = Template(SAMPLE_TEMPLATE)
        expected = [
            '    Query                                           S-TM---NGRSFRRRFRRRFTW-T-S-FRRRFRRRFS-',
            '    Peaks                                           S         S              S          S',
            '    IndividualPeaks                                   T                 T  T   ',
            '    IndividualTroughs                                  M                       ',
            '    AlphaHelix                                                                 FRRRFRRRF',
            '    AlphaHelix, AminoAcids                                     FRRRFRRRF W     ',
            '    EukaryoticLinearMotif, AminoAcids                      NGR           W     ',
        ]
        self.assertEqual(expected, match.query.noPipes)

    def testSequence(self):
        """
        The query sequence must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual('S-TM---NGRSFRRRFRRRFTW-T-S-FRRRFRRRFS-',
                         match.query.read.sequence)

    def testLandmarks(self):
        """
        The query landmarks must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(set(['AlphaHelix', 'EukaryoticLinearMotif']),
                         match.query.landmarks)

    def testTrigPoints(self):
        """
        The query trig points must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(set(['Peaks', 'IndividualPeaks',
                              'IndividualTroughs', 'AminoAcids']),
                         match.query.trigPoints)

    def testUnknownUnpairedFeatureName(self):
        """
        If a query line mentions an unknown unpaired feature name, a ValueError
        must be raised.
        """
        lines = [
            'Subject FRRRF|RRR|FRRRFW',
            'Query   STM|STW|TS',
            'Unknown  xx|   |',
        ]
        error = ("^Unknown feature name 'Unknown' found in query template$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testFeatureSequenceMismatch(self):
        """
        If a query feature sequence does not match the overall subject
        sequence, a ValueError must be raised.
        """
        lines = [
            'Subject     FRRRF|RRR|FRRRFW',
            'Query       STM|STW|TS',
            'AlphaHelix   FF|   |',
        ]
        error = ("^AlphaHelix feature sequence 'FF' found in query template "
                 "\(offset 1, length 2\) does not match the full sequence for "
                 "the query at those offsets$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)


class TestSubject(TestQueryOrSubjectMixin, TestCase):
    """
    Tests for the Subject class.
    """
    CLASS = Subject

    def testAttributes(self):
        """
        The class must have the expected attributes.
        """
        self.assertEqual('subject', self.CLASS.TYPE)
        self.assertEqual('query', self.CLASS.OTHER_TYPE)

    def testIndentLength(self):
        """
        The correct indent length must be found for the subject.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(39, match.subject.indentLength)

    def testRelevantLines(self):
        """
        The relevant lines from the template must be found for the subject.
        """
        match = Template(SAMPLE_TEMPLATE)
        expected = [
            '    Subject                            FRRRFRRRF-W------|FRRRFRRRFW---RIKR-W-W--W---|--',
            '    AlphaHelix                         FRRRFRRRF        |FRRRFRRRF                  |',
            '    AminoAcids                                   W      |         W        W W      |',
            '    EukaryoticLinearMotif, AminoAcids                   |             RIKR      W   |',
            '    AlphaHelix, EukaryoticLinearMotif                   |FRRRFRRRF    RIKR          |',
        ]
        self.assertEqual(expected, match.subject.relevantLines)

    def testFindPipes(self):
        """
        The offsets of the pipe characters must be calculated correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(56, match.subject.pipe1)
        self.assertEqual(84, match.subject.pipe2)

    def testFindPipesMissing(self):
        """
        If a subject line doesn't have 2 pipe characters, a ValueError must be
        raised.
        """
        lines = [
            'Subject FRRRF|RRR',
            'Query STM|STWTS',
        ]
        error = ("^Template line 'Subject FRRRF|RRR' in subject does not "
                 "contain two pipe characters$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testFindPipesNotAligned(self):
        """
        If the pipe characters in a subject line don't all align, a ValueError
        must be raised.
        """
        lines = [
            'Subject    FRR|RF|RRR',
            'AminoAcids   W| |',
            'Query STM|ST|WTS',
        ]
        error = ("^Template line 'AminoAcids   W| |' in subject has pipe "
                 "characters at offsets 14 and 16, but these do not match the "
                 "offsets \(14 and 17\) of the pipe characters found on the "
                 "first subject template line$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testNoPipes(self):
        """
        The relevant lines from the template must be found for the subject and
        their pipe characters must be removed.
        """
        match = Template(SAMPLE_TEMPLATE)
        expected = [
            '    Subject                            FRRRFRRRF-W------FRRRFRRRFW---RIKR-W-W--W-----',
            '    AlphaHelix                         FRRRFRRRF        FRRRFRRRF                  ',
            '    AminoAcids                                   W               W        W W      ',
            '    EukaryoticLinearMotif, AminoAcids                                RIKR      W   ',
            '    AlphaHelix, EukaryoticLinearMotif                   FRRRFRRRF    RIKR          ',
        ]
        self.assertEqual(expected, match.subject.noPipes)

    def testSequence(self):
        """
        The subject sequence must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual('FRRRFRRRF-W------FRRRFRRRFW---RIKR-W-W--W-----',
                         match.subject.read.sequence)

    def testLandmarks(self):
        """
        The subject landmarks must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(
            set(['AlphaHelix', 'EukaryoticLinearMotif']),
            match.subject.landmarks)

    def testTrigPoints(self):
        """
        The subject trig points must be extracted correctly.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(set(['AminoAcids']), match.subject.trigPoints)

    def testUnknownUnpairedFeatureName(self):
        """
        If a subject line mentions an unknown unpaired feature name, a
        ValueError must be raised.
        """
        lines = [
            'Subject FRRRF|RRR|FRRRFW',
            'Unknown    xx|   |',
            'Query   STM|STW|TS',
        ]
        error = ("^Unknown feature name 'Unknown' found in subject template$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)

    def testFeatureSequenceMismatch(self):
        """
        If a subject feature sequence does not match the overall subject
        sequence, a ValueError must be raised.
        """
        lines = [
            'Subject     FRRRF|RRR|FRRRFW',
            'AlphaHelix     FF|   |',
            'Query       STM|STW|TS',
        ]
        error = ("^AlphaHelix feature sequence 'FF' found in subject template "
                 "\(offset 3, length 2\) does not match the full sequence for "
                 "the subject at those offsets$")
        six.assertRaisesRegex(self, ValueError, error, self.CLASS, lines)


class TestBinScoreTemplate(TestCase):
    """
    Tests for the binScoreTemplate.py (in this directory) Template class.
    """

    def testUnequalPairedFeatureCount(self):
        """
        If the query and subject do not have the same number of paired
        features, Template must raise a ValueError.
        """
        template = """
            Subject                            ----|FRRRFRRRFW---RIKR-W-W--W-|
            EukaryoticLinearMotif, AminoAcids      |             RIKR      W |

            Query                              S-TM---|NGRSFRRRFRRRFTW-T-S-|
        """

        error = ('^The query and subject do not have the same number of '
                 'paired features \(0 != 1\)$')
        six.assertRaisesRegex(self, ValueError, error, Template, template)

    def testLandmarks(self):
        """
        The Template class must have the correct set of landmark names.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(set(['AlphaHelix', 'EukaryoticLinearMotif']),
                         match.landmarks)

    def testTrigPoints(self):
        """
        The Template class must have the correct set of trig point names.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(set(['Peaks', 'IndividualTroughs', 'IndividualPeaks',
                              'AminoAcids']),
                         match.trigPoints)

    def testHistogramWithOneBin(self):
        """
        The Template class must build a histogram with one bin.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(1, match.histogram.nBins)

    def testHistogram(self):
        """
        The Template class histogram must have a bin containing the correct
        query/subject pairs.
        """
        match = Template(SAMPLE_TEMPLATE)
        self.assertEqual(
            [
                {
                    'queryLandmark': Landmark(AlphaHelix.NAME,
                                              AlphaHelix.SYMBOL,
                                              offset=11, length=9),
                    'queryTrigPoint': TrigPoint(AminoAcids.NAME,
                                                AminoAcids.SYMBOL,
                                                offset=21),
                    'subjectLandmark': Landmark(EukaryoticLinearMotif.NAME,
                                                EukaryoticLinearMotif.SYMBOL,
                                                offset=30, length=4),
                    'subjectTrigPoint': TrigPoint(AminoAcids.NAME,
                                                  AminoAcids.SYMBOL,
                                                  offset=40),
                },
                {
                    'queryLandmark': Landmark(EukaryoticLinearMotif.NAME,
                                              EukaryoticLinearMotif.SYMBOL,
                                              offset=7, length=3),
                    'queryTrigPoint': TrigPoint(AminoAcids.NAME,
                                                AminoAcids.SYMBOL, offset=21),
                    'subjectLandmark': Landmark(AlphaHelix.NAME,
                                                AlphaHelix.SYMBOL,
                                                offset=17, length=9),
                    'subjectTrigPoint': Landmark(EukaryoticLinearMotif.NAME,
                                                 EukaryoticLinearMotif.SYMBOL,
                                                 offset=30, length=4)
                },
            ],
            match.histogram[0])

    def testPassDatabaseParamsThatMissLandmarksInTheTemplate(self):
        """
        If a DatabaseParameters is passed and it doesn't include all the
        landmark finders mentioned in the template, the Template
        calculateScore method must raise a ValueError.
        """
        dbParams = DatabaseParameters(landmarks=[], trigPoints=[])
        error = ('^The template mentions landmark finders \(AlphaHelix, '
                 'EukaryoticLinearMotif\) that are not present in the passed '
                 'DatabaseParameters instance$')
        match = Template(SAMPLE_TEMPLATE)
        six.assertRaisesRegex(self, ValueError, error, match.calculateScore,
                              dbParams=dbParams)

    def testPassDatabaseParamsThatMissTrigPointsInTheTemplate(self):
        """
        If a DatabaseParameters is passed and it doesn't include all the trig
        point finders mentioned in the template, the Template calculateScore
        method  must raise a ValueError.
        """
        dbParams = DatabaseParameters(
            landmarks=['AlphaHelix', 'EukaryoticLinearMotif'], trigPoints=[])
        error = ('^The template mentions trig point finders \(AminoAcids, '
                 'IndividualPeaks, IndividualTroughs, Peaks\) that are not '
                 'present in the passed DatabaseParameters instance$')
        match = Template(SAMPLE_TEMPLATE)
        six.assertRaisesRegex(self, ValueError, error, match.calculateScore,
                              dbParams=dbParams)

    def testPassUnknownBinScoreMethod(self):
        """
        If a FindParameters is passed and its bin score method is unknown,
        the Template calculateScore method must raise a ValueError.
        """
        findParams = FindParameters(binScoreMethod='unknown')
        error = ("^Unknown bin score method 'unknown'$")
        match = Template(SAMPLE_TEMPLATE)
        six.assertRaisesRegex(self, ValueError, error, match.calculateScore,
                              findParams=findParams)
