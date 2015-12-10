import six
from unittest import TestCase

from light.utils import maxWithDefault, minWithDefault


class TestMaxWithDefault(TestCase):
    """
    Tests for light.utils.maxWithDefault function.
    """

    def testMaxWithEmptyListAndNoDefault(self):
        """
        If maxWithDefault is called with an empty list and no default, None
        must be returned.
        """
        error = '^max\(\) arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error, maxWithDefault, [])

    def testMaxWithEmptyListAndDefault(self):
        """
        If maxWithDefault is called with an empty list and a default value,
        then the default value must be returned.
        """
        self.assertIs(3, maxWithDefault([], default=3))

    def testMaxWithNonEmptyListAndNoDefault(self):
        """
        If maxWithDefault is called with a non-empty list (and no default),
        then the max of the list must be returned.
        """
        self.assertIs(4, maxWithDefault([3, 4, 2]))

    def testMaxWithNonEmptyListAndDefault(self):
        """
        If maxWithDefault is called with a non-empty list (and a default),
        then the max of the list must be returned.
        """
        self.assertIs(4, maxWithDefault([3, 4, 2], default=7))


class TestMinWithDefault(TestCase):
    """
    Tests for light.utils.minWithDefault function.
    """

    def testMinWithEmptyListAndNoDefault(self):
        """
        If minWithDefault is called with an empty list and no default, a
        ValueError must be raised.
        """
        error = '^min\(\) arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error, minWithDefault, [])

    def testMinWithEmptyListAndDefault(self):
        """
        If minWithDefault is called with an empty list and a default value,
        then the default value must be returned.
        """
        self.assertIs(3, minWithDefault([], default=3))

    def testMinWithNonEmptyListAndNoDefault(self):
        """
        If minWithDefault is called with a non-empty list (and no default),
        then the min of the list must be returned.
        """
        self.assertIs(2, minWithDefault([3, 4, 2]))

    def testMinWithNonEmptyListAndDefault(self):
        """
        If minWithDefault is called with a non-empty list (and a default),
        then the min of the list must be returned.
        """
        self.assertIs(2, minWithDefault([3, 4, 2], default=7))
