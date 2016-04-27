import six
from unittest import TestCase

from light.performance.utils import pdbNameToPythonName


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
