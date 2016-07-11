import six
from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from ..mocking import mockOpen

from light.performance.pdb import (
    pdbNameToPythonName, pythonNameToPdbName, loadObsoletePDB)


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


class TestLoadObsoletePDB(TestCase):
    """
    Tests for the loadObsoletePDB function.
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
                                  loadObsoletePDB, 'filename')

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
                                  loadObsoletePDB, 'filename')

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
                                  loadObsoletePDB, 'filename')

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
                                  loadObsoletePDB, 'filename')

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
            result = loadObsoletePDB('filename')
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
            result = loadObsoletePDB('filename')
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
