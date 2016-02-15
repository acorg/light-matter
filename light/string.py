class MultilineString:
    """
    Incrementally build a multi-line indented string.

    @param indent: A C{str} indent that will be inserted (zero or more times)
        at the start of each line, depending on the current indent level
        multiplier.
    @param margin: A C{str} that will always be inserted once at the start of
        each line.
    """
    def __init__(self, indent='  ', margin=''):
        self._indent = indent
        self._margin = margin
        self._level = 0
        self._strings = []

    def indent(self):
        """
        Increase the number of times the indent is inserted before each line.

        @return: self, to allow our caller to chain calls.
        """
        self._level += 1
        return self

    def outdent(self):
        """
        Decrease the number of times the indent is inserted before each line.

        @return: self, to allow our caller to chain calls.
        """
        self._level -= 1
        if self._level < 0:
            self._level = 0
        return self

    def append(self, s, verbatim=False):
        """
        Add a line to the multi-line string.

        @param s: A C{str} to add to the multi-line string.
        @param verbatim: If C{True}, C{s} will be split on '\n' and the
            resulting strings will be added to this multi-line string as-is
            (with no margin or indent). This is useful for adding a
            pre-formatted multi-line string (perhaps obtained via a network
            call).
        @return: self, to allow our caller to chain calls.
        """
        if verbatim:
            self._strings.extend(s.split('\n'))
        else:
            self._strings.append(self._margin + self._level * self._indent + s)
        return self

    def extend(self, strings):
        """
        Add a multiple lines to the multi-line string.

        @param strings: An iterable of C{str} instances to add to the
            multi-line string.
        @return: self, to allow our caller to chain calls.
        """
        indent = self._margin + self._level * self._indent
        self._strings.extend(indent + string for string in strings)
        return self

    def __str__(self):
        """
        Get a version of the string suitable for printing.

        @return: A single C{str} of the multi-line string with embedded
            newlines (but no final newline).
        """
        return '\n'.join(self._strings)

    def __len__(self):
        """
        Get the length of a multi-line string. Note that this can be used
        for truth testing as well (using 'if' and 'if not').

        @return: An C{int} length of the string.
        """
        # Note that with the following simplistic implementation it is
        # inefficient to write e.g., 'return str(s) if s else None' when s
        # is a MultilineString as that converts s to a str twice. It would
        # be more efficient to use 'x = str(s); return x if x else None'.
        # To alleviate this we could instead provide a method that tested
        # the length of self._strings. But for now, be dumb.
        return len(str(self))

    def lineCount(self):
        """
        How many lines are in this multi-line string?

        @return: The C{int} number of lines in this multi-line string.
        """
        return len(self._strings)
