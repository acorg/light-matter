from os import mkdir
from os.path import exists, isdir, join


def makeDir(path):
    """
    Check to see if a path exists and whether it's a directory. Create it (as
    a dir) if it's non-existent.

    @raise RuntimeError: if C{path} exists but is not a directory.
    @param path: The C{str} path for the directory.
    """
    if exists(path):
        if not isdir(path):
            raise RuntimeError('Output path %r already exists but is not a '
                               'directory.' % path)
    else:
        mkdir(path)


def makeOutputDir(*path):
    """
    Create an output directory if it doesn't already exist.

    @param path: A C{list} of C{str} path components to add to
        C{testArgs.outputDir}. The full path to the directory is made, as in
        mkdir -p at the shell.

    @return: The C{str} path to the output directory.
    """
    # This import cannot be done at the top level, it needs to be deferred
    # until light/performance/bin/perf.py sets it.
    from light.performance import testArgs
    result = testArgs.outputDir

    for p in path:
        result = join(result, p)
        makeDir(result)

    return result
