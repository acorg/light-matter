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


def pdbNameToPythonName(name, raiseOnError=True):
    """
    Convert a PDB protein name to a name that can be used in Python
    (and the filesystem).

    @param: A C{str} PDB name, usually in the form "NXXX:X" where N is a digit,
        and each X is an ASCII letter. Other forms are possible, for example
        with a trailing ":sequence" or "NXXX-X".
    @param raiseOnError: When a name cannot be converted: if C{True} raise a
        C{RuntimeError}, else return the original name.
    @raise RuntimeError: If the format of C{name} isn't recognized and
        C{raiseOnError} is C{True}.
    @return A lowercase C{str} name in the form "pdb_NXXX_X" (to follow the
        above example) that can be used as a Python variable, a name that can
        be used for a file to be imported as a Python module.
    """
    origName = name
    name = name.lower()

    if name.endswith(':sequence'):
        name = name[:-9]

    n = len(name)

    if n == 6 and name[4] in (':', '-', '_'):
        pdbId = name[0:4]
        suffix = name[5]
    else:
        if raiseOnError:
            raise RuntimeError(
                'Could not convert PDB name %r to a Python name' % origName)
        else:
            return origName

    return 'pdb_' + pdbId + '_' + suffix
