import six
from unittest import TestCase

from light.performance.utils import pdbNameToPythonName, pythonNameToPdbName


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
