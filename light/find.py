from light.statistics import HasAllBases

ALLSTATS = [HasAllBases]
DNASTATS = []
RNASTATS = []


def find(which):
    """
    A function to find all subclasses that should be executed. To be used to
    find the statistics to be run.

    TODO: Needs to be updated with more specific filtering methods.
          Expand globals!

    @param which: The type of class that should be found. Can eigher be 'all',
        'dna', 'rna'.
    """
    if which == 'all':
        return ALLSTATS
    elif which == 'dna':
        return DNASTATS
    elif which == 'rna':
        return RNASTATS
    else:
        return False
