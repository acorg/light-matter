import six
from unittest import TestCase

from light.performance.utils import (
    convertDictKeysFromPDBToPython, pdbNameToPythonName)


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


class TestConvertDictKeysFromPDBToPython(TestCase):
    """
    Tests for the light.performance.utils.convertDictKeysFromPDBToPython
    function.
    """
    def testEmpty(self):
        """
        Passing an empty dict must result in an empty dict.
        """
        self.assertEqual({}, convertDictKeysFromPDBToPython({}))

    def testUnrecognizable(self):
        """
        An unrecognizable name in a dict must result in a RuntimeError.
        """
        error = "^Could not convert PDB name 'xxxxxxx' to a Python name$"
        six.assertRaisesRegex(self, RuntimeError, error,
                              convertDictKeysFromPDBToPython,
                              {'xxxxxxx': 44})

    def testUnrecognizableWithoutRaising(self):
        """
        An unrecognizable name in a dict must result in the original name
        being in the result (not in the raising of a RuntimeError) if
        raiseOnError is passed as False.
        """
        d = {'x' * 10: 33}
        self.assertEqual(d, convertDictKeysFromPDBToPython(
            d, raiseOnError=False))

    def testCanonical(self):
        """
        Passing a dict with canonical PDB names must return as expected.
        """
        self.assertEqual(
            dict.fromkeys(['pdb_1mpl_a', 'pdb_2abc_d']),
            convertDictKeysFromPDBToPython(
                dict.fromkeys(['1MPL:A', '2ABC:D'])))

    def testValuesArePreserverd(self):
        """
        The values in the result must be those in the argument.
        """
        o = object()
        result = convertDictKeysFromPDBToPython({'1MPL:A': o})
        self.assertIs(o, result['pdb_1mpl_a'])
