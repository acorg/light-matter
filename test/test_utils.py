from unittest import TestCase

from light.utils import convertAAToAAProperties, HYDROPHOBIC, HYDROPHILIC


class TestConvertAAToAAProperties(TestCase):
    """
    Test for the utils.convertAAToAAProperties function.
    """
    def testAssertNumberOfProperties(self):
        """
        Test that an assertion is raised if not exactly two properties are
        specified.
        """
        sequence = 'ASDGEBHSDTDSCV'
        self.assertRaisesRegexp(Exception, 'Two properties must be specified',
                                convertAAToAAProperties, sequence, [1, 2, 3])

    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        sequence = 'ASDGEBHSDTDSCV'
        result = convertAAToAAProperties(sequence,
                                         [HYDROPHOBIC, HYDROPHILIC])
        self.assertEqual('IOOIOIOOIOOII', result)
