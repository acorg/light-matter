from functools import total_ordering

from light.string import MultilineString


@total_ordering
class Finder(object):
    """
    Superclass of all light.landmark and light.trig feature finders.

    @param dbParams: A C{DatabaseParameters} instance.
    """

    # PARAMETERS can be used by Finder subclasses to define finder-specific
    # parameters. Its keys are parameter names, and its values are
    # dictionaries suitable for passing to the add_argument method of an
    # argparse.ArgumentParser instance (see addArgsToParser below). See
    # landmarks/random.py for an example usage.
    PARAMETERS = {}

    def __eq__(self, other):
        """
        Should two finders compare equal?

        Note that no distinction is made between landmark and trig point
        finders.

        @param: Another C{Finder} instance.
        @return: C{True} if the two finders are the same.
        """
        return (self.NAME, self.SYMBOL) == (other.NAME, other.SYMBOL)

    def __lt__(self, other):
        """
        Do we compare less than another finder?

        Note that no distinction is made between landmark and trig point
        finders.

        @param: Another C{Finder} instance.
        @return: C{True} if C{self} is less than the other finder.
        """
        return (self.NAME, self.SYMBOL) < (other.NAME, other.SYMBOL)

    def addArgsToParser(self, parser):
        """
        Add this finder's parameter arguments (if any) to an argparse parser.

        @param parser: An C{argparse.ArgumentParser} instance.
        """
        for parameterName, parameterKwds in self.PARAMETERS.items():
            parser.add_argument('--' + parameterName, **parameterKwds)

    def state(self):
        """
        Get the state of the finder (its parameter settings).

        @return: A C{dict} of parameter name, parameter value.
        """
        result = {}
        for parameterName in self.PARAMETERS:
            result[parameterName] = getattr(self, parameterName)
        return result

    def initializeFromArgs(self, args):
        """
        Set a finder's parameters from values given on the command line.

        @param args: Command line arguments as returned by the C{argparse}
            C{parse_args} method.
        """
        for parameterName, parameterKwds in self.PARAMETERS.items():
            setattr(self, parameterName,
                    getattr(args, parameterName, parameterKwds['default']))

    def initializeFromKeywords(self, **kwargs):
        """
        Set all a finder's parameters from keyword arguments (if present in
        C{kwargs}) or using the finder's default value (if not).

        @param kwargs: A C{dict} of keywords and values.
        """
        print('  init', self.NAME)
        for parameterName, parameterKwds in self.PARAMETERS.items():
            print('      ', parameterName, parameterKwds)
            setattr(self, parameterName,
                    kwargs.get(parameterName, parameterKwds['default']))

    def print_(self, margin=''):
        """
        Create a C{str} with the finder name and its parameter values (if any).

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the finder and its parameters.
        """
        result = MultilineString(margin=margin)
        result.append(self.NAME)
        if self.PARAMETERS:
            result.indent()
            result.extend([
                '%s: %s' % (parameterName, getattr(self, parameterName))
                for parameterName in sorted(self.PARAMETERS)])
            result.outdent()
        return str(result)
