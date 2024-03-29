import six
from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from ..mocking import mockOpen

from light.performance.pdb import (
    pdbNameToPythonName, pythonNameToPdbName, loadObsolete, loadResolution,
    loadEntries)


class TestPdbNameToPythonName(TestCase):
    """
    Tests for the light.performance.utils.pdbNameToPythonName function.
    """
    def testUnrecognizable(self):
        """
        An unrecognizable name must result in a RuntimeError.
        """
        error = "^Could not convert PDB name 'xxxxxxx' to a Python name$"
        six.assertRaisesRegex(self, RuntimeError, error,
                              pdbNameToPythonName, 'xxxxxxx')

    def testUnrecognizableWithoutRaising(self):
        """
        An unrecognizable name must return the original name (not raise
        a RuntimeError) if raiseOnError is passed as False.
        """
        self.assertEqual('x' * 10,
                         pdbNameToPythonName('x' * 10, raiseOnError=False))

    def testCanonicalPDBUpperCase(self):
        """
        A canonical DXXX:X PDB name must be converted as expected.
        """
        self.assertEqual('pdb_1mla_a', pdbNameToPythonName('1MLA:A'))

    def testCanonicalPDBMixedCase(self):
        """
        A canonical DXXX:X PDB name in mixed case must be converted as
        expected.
        """
        self.assertEqual('pdb_1mla_a', pdbNameToPythonName('1MLa:a'))

    def testCanonicalPDBLowerCase(self):
        """
        A canonical DXXX:X PDB name in lower case must be converted as
        expected.
        """
        self.assertEqual('pdb_1mla_a', pdbNameToPythonName('1mla:a'))

    def testCanonicalPDBWithHyphen(self):
        """
        A canonical DXXX-X PDB name must be converted as expected.
        """
        self.assertEqual('pdb_1mla_a', pdbNameToPythonName('1MLA-A'))

    def testCanonicalPDBWithUnderscore(self):
        """
        A canonical DXXX_X PDB name must be converted as expected.
        """
        self.assertEqual('pdb_1mla_a', pdbNameToPythonName('1MLA_A'))

    def testCanonicalPDBWithColonSequenceSuffix(self):
        """
        A DXXX:X:sequence PDB name must be converted as expected.
        """
        self.assertEqual('pdb_1mla_a',
                         pdbNameToPythonName('1MLA:A:sequence'))


class TestPythonNameToPdbName(TestCase):
    """
    Tests for the light.performance.utils.pythonNameToPdbName function.
    """
    def testUnrecognizable(self):
        """
        Passing an unrecognizable name must result in the original name.
        """
        self.assertEqual('x' * 10, pythonNameToPdbName('x' * 10))

    def testUnrecognizableDueToSuffix(self):
        """
        Passing a name that has a valid prefix but which is too long must
        result in the original name.
        """
        self.assertEqual('pdb_1mla_aa', pythonNameToPdbName('pdb_1mla_aa'))

    def testUnrecognizableDueToPreix(self):
        """
        Passing a name that has a valid suffix but which is too long must
        result in the original name.
        """
        self.assertEqual('apdb_1mla_a', pythonNameToPdbName('apdb_1mla_a'))

    def testRecognizedLowerCase(self):
        """
        Passing an recognizable Python name must result in the expected PDB
        name.
        """
        self.assertEqual('1MLA:A', pythonNameToPdbName('pdb_1mla_a'))


class TestLoadObsolete(TestCase):
    """
    Tests for the loadObsolete function.
    """

    FIRST_LINE = ' LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS'

    def testIncorrectFirstLine(self):
        """
        A ValueError must be raised if the first line is not as expected.
        """
        mockOpener = mockOpen(read_data='xxx\n')
        with patch.object(builtins, 'open', mockOpener):
            error = "^Unexpected first line of input file 'filename'."
            six.assertRaisesRegex(self, ValueError, error,
                                  loadObsolete, 'filename')

    def testMissingOBSLTE(self):
        """
        A ValueError must be raised if an input line does not start with
        OBSLTE.
        """
        data = '\n'.join([
            self.FIRST_LINE,
            'xxx 1 2',
        ]) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Line 2 of input file 'filename' does not start with "
                     "OBSLTE\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadObsolete, 'filename')

    def testTooFewFields(self):
        """
        A ValueError must be raised if an input line has too few fields.
        """
        data = '\n'.join([
            self.FIRST_LINE,
            'OBSLTE 1',
        ]) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Line 2 of input file 'filename' has 2 fields "
                     "\(expected 3 or more\)\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadObsolete, 'filename')

    def testRepeatedId(self):
        """
        A ValueError must be raised if an input line repeats a PDB id.
        """
        data = '\n'.join([
            self.FIRST_LINE,
            'OBSLTE date 4HLA',
            'OBSLTE date 4HLA',
        ]) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Repeated PDB id '4HLA' found on line 3 of input "
                     "file 'filename'\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadObsolete, 'filename')

    def testNoSuccessors(self):
        """
        If no successors are given for a sequence id, the successors list
        in the result must be empty.
        """
        data = '\n'.join([
            self.FIRST_LINE,
            'OBSLTE date 4HLA',
        ]) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadObsolete('filename')
        self.assertEqual([], result['4HLA']['successors'])

    def testExpectedCorrectResult(self):
        """
        If a valid input file is passed, the expected result must be returned.
        """
        data = '\n'.join([
            self.FIRST_LINE,
            'OBSLTE date1 4HLA 5HML',
            'OBSLTE date2 2ABC 3DEF 4GHI',
        ]) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadObsolete('filename')
        self.assertEqual(
            {
                '2ABC': {
                    'date': 'date2',
                    'id': '2ABC',
                    'successors': ['3DEF', '4GHI'],
                },
                '4HLA': {
                    'date': 'date1',
                    'id': '4HLA',
                    'successors': ['5HML']
                }
            },
            result)


class TestLoadResolution(TestCase):
    """
    Tests for the loadResolution function.
    """

    HEADER_LINES = [
        'PROTEIN DATA BANK LIST OF IDCODES AND DATA RESOLUTION VALUES',
        'Fri Jul 01 14:17:44 EDT 2016',  # Date line - not checked.
        ('RESOLUTION VALUE IS -1.00 FOR ENTRIES DERIVED FROM NMR AND OTHER '
         'EXPERIMENT METHODS (NOT INCLUDING X-RAY) IN WHICH THE FIELD '
         'REFINE.LS_D_RES_HIGH IS EMPTY'),
        '',
        'IDCODE       RESOLUTION',
        '------  -    ----------',
    ]

    def testIncorrectLine1(self):
        """
        A ValueError must be raised if the first line is not as expected.
        """
        mockOpener = mockOpen(read_data='xxx\n')
        with patch.object(builtins, 'open', mockOpener):
            error = "^Line 1 of 'filename' was expected to be 'PROTEIN DATA "
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testIncorrectLine3(self):
        """
        A ValueError must be raised if the third line is not as expected.
        """
        data = '\n'.join(self.HEADER_LINES[:2] + ['xxx']) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = (
                "^Line 3 of 'filename' was expected to be 'RESOLUTION VALUE ")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testIncorrectLine4(self):
        """
        A ValueError must be raised if the fourth line is not as expected.
        """
        data = '\n'.join(self.HEADER_LINES[:3] + ['xxx']) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^Line 4 of 'filename' was expected to be '"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testIncorrectLine5(self):
        """
        A ValueError must be raised if the fifth line is not as expected.
        """
        data = '\n'.join(self.HEADER_LINES[:4] + ['xxx']) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^Line 5 of 'filename' was expected to be 'IDCODE"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testIncorrectLine6(self):
        """
        A ValueError must be raised if the sixth line is not as expected.
        """
        data = '\n'.join(self.HEADER_LINES[:5] + ['xxx']) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^Line 6 of 'filename' was expected to be '------"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testMissingSemicolon(self):
        """
        A ValueError must be raised if an input line does not have a semicolon
        separator.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id + 5.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Line 7 of 'filename' does not contain expected "
                     "semicolon separator\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testTooFewFields(self):
        """
        A ValueError must be raised if an input line has too few fields.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id 3.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = '^not enough values to unpack \(expected 3, got 2\)$'
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename')

    def testRepeatedIdBest(self):
        """
        The best (lowest) resolution must be returned if an input line repeats
        a PDB id and the conflict resolution is set to 'best'.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id1 ; 3.0', 'id1 ; 2.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadResolution('filename', 'best')
        self.assertEqual(
            {
                'id1': 2.0,
            },
            result)

    def testRepeatedIdWorst(self):
        """
        The best (highest) resolution must be returned if an input line repeats
        a PDB id and the conflict resolution is set to 'worst'.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id1 ; 3.0', 'id1 ; 2.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadResolution('filename', 'worst')
        self.assertEqual(
            {
                'id1': 3.0,
            },
            result)

    def testRepeatedIdRaise(self):
        """
        A ValueError must be raised if an input line repeats a PDB id and
        the conflict resolution is set to 'raise'.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id1 ; 3.0', 'id1 ; 3.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Repeated PDB id 'id1' found on line 8 of input "
                     "file 'filename'\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename', 'raise')

    def testRepeatedIdUnknownResolutionMethod(self):
        """
        A ValueError must be raised if an input line repeats a PDB id and
        the conflict resolution is unknown.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id1 ; 3.0', 'id1 ; 3.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^whenConflicting must be one of best, raise, worst\.$"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadResolution, 'filename', 'unknown')

    def testExpectedCorrectResult(self):
        """
        If a valid input file is passed, the expected result must be returned.
        """
        data = '\n'.join(
            self.HEADER_LINES + ['id1 ; 3.0', 'id2 ; 2.0', 'id3 ; -1.0']
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadResolution('filename')
        self.assertEqual(
            {
                'id1': 3.0,
                'id2': 2.0,
                'id3': -1.0,
            },
            result)


class TestLoadEntries(TestCase):
    """
    Tests for the loadEntries function.
    """

    HEADER_LINES = [
        'IDCODE, HEADER, ACCESSION DATE, COMPOUND, SOURCE, AUTHOR LIST, '
        'RESOLUTION, EXPERIMENT TYPE (IF NOT X-RAY)',
        '------- ------- --------------- --------- ------- ------------ '
        '----------- ---------------------------------------------------'
        '---------------------------------------------------------------'
        '---------------------------------------------------------------'
        '---------------------------------------------------------------'
        '----------------------------------------------------------'
    ]

    def testIncorrectLine1(self):
        """
        A ValueError must be raised if the first line is not as expected.
        """
        mockOpener = mockOpen(read_data='xxx\n')
        with patch.object(builtins, 'open', mockOpener):
            error = "^Line 1 of 'filename' was expected to be 'IDCODE, "
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testIncorrectLine2(self):
        """
        A ValueError must be raised if the second line is not as expected.
        """
        data = '\n'.join(self.HEADER_LINES[:1] + ['xxx']) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = (
                "^Line 2 of 'filename' was expected to be '------- -------")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testTooFewFields(self):
        """
        A ValueError must be raised if an input line has too few fields.
        """
        data = '\n'.join(
            self.HEADER_LINES + ["1\t2\t3"]
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = '^not enough values to unpack \(expected 8, got 3\)$'
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testRepeatedId(self):
        """
        A ValueError must be raised if an input line repeats a PDB id and
        the conflict resolution is set to 'raise'.
        """
        data = '\n'.join(
            self.HEADER_LINES + [
                "100D\tDNA\t12/05/94\tcompound\tsource\tauths\tres\ttype",
                "100D\tDNA\t01/27/03\tcompound\tsource\tauths\tres\ttype",
            ]
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("^Repeated PDB id '100D' found on line 4 of input "
                     "file 'filename'\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testUnparseableDate(self):
        """
        A ValueError must be raised if an input line contains a date that
        cannot be split by '/'.
        """
        data = '\n'.join(
            self.HEADER_LINES + [
                "100D\tDNA\t12-05-94\tcompound\tsource\tauths\tres\ttype",
            ]
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^invalid literal for int\(\) with base 10: '12-05-94'$"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testNotEnoughDateFields(self):
        """
        A ValueError must be raised if an input line contains a date that
        cannot be split into 3 parts by '/'.
        """
        data = '\n'.join(
            self.HEADER_LINES + [
                "100D\tDNA\t12/05\tcompound\tsource\tauths\tres\ttype",
            ]
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^not enough values to unpack \(expected 3, got 2\)$"
            six.assertRaisesRegex(self, ValueError, error,
                                  loadEntries, 'filename')

    def testExpectedCorrectResult(self):
        """
        If a valid input file is passed, the expected result must be returned.
        """
        data = '\n'.join(
            self.HEADER_LINES + [
                "100D\tDNA\t12/05/94\tcompound\tsource\tauths\tres\ttype",
                "103F\tDNA\t01/27/03\tcompound\tsource\tauths\tres\ttype",
                "123G\tDNA\t1/1/00\tcompound\tsource\tauths\tres\ttype",
            ]
        ) + '\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            result = loadEntries('filename')
        self.assertEqual(
            {
                '100D': {
                    'day': 5,
                    'month': 12,
                    'year': 1994,
                },
                '103F': {
                    'day': 27,
                    'month': 1,
                    'year': 2003,
                },
                '123G': {
                    'day': 1,
                    'month': 1,
                    'year': 2000,
                },
            },
            result)
