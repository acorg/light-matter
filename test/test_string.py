from unittest import TestCase

from light.string import MultilineString


class TestMultilineString(TestCase):
    """
    Tests for the light.string.MultilineString class.
    """
    def testEmpty(self):
        """
        If a multiline string has nothing added to it, its str representation
        must be the empty string.
        """
        self.assertEqual('', str(MultilineString()))

    def testEmptyIsFalse(self):
        """
        If a multiline string has nothing added to it, it must be considered
        False in a Boolean test.
        """
        self.assertFalse(MultilineString())

    def testLenZero(self):
        """
        If a multiline string has nothing added to it, its length must be zero.
        """
        self.assertEqual(0, len(MultilineString()))

    def testAppendOneLine(self):
        """
        If a multiline string has one thing appended to it, its str
        representation must be just that string.
        """
        s = MultilineString()
        s.append('test')
        self.assertEqual('test', str(s))

    def testOneLineTrue(self):
        """
        If a multiline string has one thing appended to it, it must be
        considered True in a Boolean test.
        """
        s = MultilineString()
        s.append('test')
        self.assertTrue(s)

    def testOneLineLength(self):
        """
        If a multiline string has one thing appended to it, its length
        must be the length of that string.
        """
        s = MultilineString()
        s.append('test')
        self.assertEqual(4, len(s))

    def testAppendReturnsSelf(self):
        """
        The append method must return the MultilineString instance.
        """
        s = MultilineString()
        self.assertIs(s, s.append('test'))

    def testAppendTwoLines(self):
        """
        If a multiline string has two things appended to it, its str
        representation must have the two strings, separated by a newline.
        """
        s = MultilineString()
        s.append('test1')
        s.append('test2')
        self.assertEqual('test1\ntest2', str(s))

    def testTwoLinesLength(self):
        """
        If a multiline string has two things appended to it, its length
        must be the sum of the two lengths plus one for the newline separator.
        """
        s = MultilineString()
        s.append('123')
        s.append('45678')
        self.assertEqual(9, len(s))

    def testAppendTwoLinesWithMargin(self):
        """
        If a multiline string with a left margin has two things appended to
        it, its str representation must have the two strings, separated by
        a newline and each line must start with the margin.
        """
        s = MultilineString(margin='+ ')
        s.append('test1')
        s.append('test2')
        self.assertEqual('+ test1\n+ test2', str(s))

    def testExtendReturnsSelf(self):
        """
        The extend method must return the MultilineString instance.
        """
        s = MultilineString()
        self.assertIs(s, s.extend(['test']))

    def testExtendTwoLines(self):
        """
        If a multiline string is extended by two strings, its str
        representation must have the two strings, separated by a newline.
        """
        s = MultilineString()
        s.extend(['test1', 'test2'])
        self.assertEqual('test1\ntest2', str(s))

    def testIndentReturnsSelf(self):
        """
        The indent method must return the MultilineString instance.
        """
        s = MultilineString()
        self.assertIs(s, s.indent())

    def testDefaultIndent(self):
        """
        If a multiline string has two things appended to it, and the second is
        indented, its str representation must have the two strings, separated
        by a newline, and the second string must be indented by two spaces (the
        default).
        """
        s = MultilineString()
        s.append('test1')
        s.indent()
        s.append('test2')
        self.assertEqual('test1\n  test2', str(s))

    def testNonDefaultIndent(self):
        """
        If a multiline string has two things appended to it, and the second is
        indented, its str representation must have the two strings, separated
        by a newline, and the second string must be indented using the indent
        value passed to __init__.
        """
        s = MultilineString(indent='\t')
        s.append('test1')
        s.indent()
        s.append('test2')
        self.assertEqual('test1\n\ttest2', str(s))

    def testOutdentReturnsSelf(self):
        """
        The outdent method must return the MultilineString instance.
        """
        s = MultilineString()
        self.assertIs(s, s.outdent())

    def testIndentOutdent(self):
        """
        If a multiline string has three things appended to it, and only the
        second is indented, its str representation must have the three strings,
        separated by newlines, and the second string must be indented.
        """
        s = MultilineString()
        s.append('test1')
        s.indent()
        s.append('test2')
        s.outdent()
        s.append('test3')
        self.assertEqual('test1\n  test2\ntest3', str(s))

    def testIndentOutdentWithMargin(self):
        """
        If a multiline string with a margin has three things appended to it,
        and only the second is indented, its str representation must have the
        three strings, separated by newlines, the second string must be
        indented, and each line must be preceeded by the margin.
        """
        s = MultilineString(indent='\t', margin='>>> ')
        s.append('test1')
        s.indent()
        s.append('test2')
        s.outdent()
        s.append('test3')
        self.assertEqual('>>> test1\n>>> \ttest2\n>>> test3', str(s))

    def testAppendVerbatim(self):
        """
        If a multiline string is given a verbatim string to append, its
        representation must be as expected.
        """
        s = MultilineString()
        s.indent()
        s.append('test1')
        s.append('      test2\n      test3', verbatim=True)
        s.append('test4')
        self.assertEqual('  test1\n      test2\n      test3\n  test4', str(s))
