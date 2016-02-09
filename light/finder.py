from functools import total_ordering


@total_ordering
class Finder(object):
    """
    Superclass of all light.landmark and light.trig feature finders.

    @param dbParams: A C{DatabaseParameters} instance.
    """
    def __init__(self, dbParams=None):
        if dbParams:
            self._dbParams = dbParams
        else:
            # We can't import DatabaseParameters at the top level as it
            # causes a circular import problem.
            from light.parameters import DatabaseParameters
            self._dbParams = DatabaseParameters()

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
