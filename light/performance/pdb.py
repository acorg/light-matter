def pdbNameToPythonName(name, raiseOnError=True):
    """
    Convert a PDB protein name to a name that can be used in Python
    (and the filesystem).

    See the corresponding pythonNameToPdbName function below.

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


def pythonNameToPdbName(name):
    """
    Convert a Python name to a PDB name.

    See the corresponding pdbNameToPythonName function above.

    @param name: A C{str} in the form pdb_NXXX_X (e.g., pdb_2hla_a).
    @return: A C{str} PDB name in the form NXXX:X if C{name} is recognized
        as a Python name (as above), else just return C{name}.
    """
    if (len(name) == 10 and name.startswith('pdb_') and name[8] == '_'):
        return (name[4:8] + ':' + name[9]).upper()
    else:
        return name


def loadObsolete(filename):
    """
    Read a file of obsolete PDB structures.

    The file will typically be the one at
    ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat (we currently have a
    copy of this in data/pdb-20160711-obsolete.txt).

    @param filename: A C{str} file name to read.
    @raises ValueError: If the input file is not in the expected format.
    @return: A C{dict} whose keys are C{str} PDB identifiers and whose
        values are C{dict}s with the following keys:

            id: A C{str} PDB identifier, as found in the file (i.e., no
                case conversion is done).
            date: A C{str}, currently in the format e.g., 15-OCT-95.
            successors: A C{list} of C{str} PDB ids (may be empty).
    """
    result = {}
    firstLine = ' LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS\n'
    with open(filename) as fp:
        first = True
        for lineNum, line in enumerate(fp, start=1):
            if first:
                if line != firstLine:
                    raise ValueError(
                        'Unexpected first line of input file %r. '
                        '%r does not match %r.' % (filename, line, firstLine))
                first = False
            else:
                fields = line.split()
                if len(fields) == 3:
                    obsolete, date, pdbId = fields
                    successors = []
                elif len(fields) > 3:
                    obsolete, date, pdbId = fields[0:3]
                    successors = fields[3:]
                else:
                    raise ValueError(
                        'Line %d of input file %r has %d fields (expected 3 '
                        'or more).' % (lineNum, filename, len(fields)))

                if obsolete != 'OBSLTE':
                    raise ValueError(
                        'Line %d of input file %r '
                        'does not start with OBSLTE.' % (lineNum, filename))

                if pdbId in result:
                    raise ValueError(
                        'Repeated PDB id %r found on line %d of input file %r.'
                        % (pdbId, lineNum, filename))

                result[pdbId] = {
                    'id': pdbId,
                    'date': date,
                    'successors': successors,
                }

    return result


def loadResolution(filename, whenConflicting='best'):
    """
    Read a file of PDB structure resolutions.

    The file will typically be the one at
    ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx (we currently
    have a copy of this in data/pdb-20160711-resolution.txt).

    @param filename: A C{str} file name to read.
    @param whenConflicting: A C{str} value indicating what to do when a
        structure is found more than once on input. Possible values are
        'best', 'worst', 'raise', to return the best (numerically smallest)
        resolution, worst (numerically highest) resolution, or to raise
        a C{ValueError}.
    @raises ValueError: If the input file is not in the expected format or
        if C{whenConflicting} is invalid.
    @return: A C{dict} whose keys are C{str} PDB identifiers and whose
        values are C{float} resolutions (in Angstroms).
    """
    conflictResolutions = ('best', 'worst', 'raise')
    if whenConflicting not in conflictResolutions:
        raise ValueError('whenConflicting must be one of %s.' %
                         ', '.join(sorted(conflictResolutions)))

    result = {}

    headerLines = [
        'PROTEIN DATA BANK LIST OF IDCODES AND DATA RESOLUTION VALUES\n',
        None,  # Date line - not checked.
        ('RESOLUTION VALUE IS -1.00 FOR ENTRIES DERIVED FROM NMR AND OTHER '
         'EXPERIMENT METHODS (NOT INCLUDING X-RAY) IN WHICH THE FIELD '
         'REFINE.LS_D_RES_HIGH IS EMPTY\n'),
        '\n',
        'IDCODE       RESOLUTION\n',
        '------  -    ----------\n',
    ]

    headerLineCount = len(headerLines)

    with open(filename) as fp:
        for lineNum, line in enumerate(fp, start=1):
            if lineNum <= headerLineCount:
                expected = headerLines[lineNum - 1]
                if expected is not None:
                    if line != expected:
                        raise ValueError(
                            'Line %d of %r was expected to be %r but was %r.'
                            % (lineNum, filename, expected, line))
            else:
                pdbId, semicolon, resolution = line.split()
                if semicolon != ';':
                    raise ValueError(
                        'Line %d of %r does not contain expected semicolon '
                        'separator.' % (lineNum, filename))

                resolution = float(resolution)

                try:
                    existingResolution = result[pdbId]
                except KeyError:
                    result[pdbId] = resolution
                else:
                    if whenConflicting == 'best':
                        if resolution < existingResolution:
                            result[pdbId] = resolution
                    elif whenConflicting == 'worst':
                        if resolution > existingResolution:
                            result[pdbId] = resolution
                    else:
                        raise ValueError(
                            'Repeated PDB id %r found on line %d of input '
                            'file %r.' % (pdbId, lineNum, filename))

    return result


def loadEntries(filename):
    """
    Read a file of PDB entries information.

    The file will typically be the one at
    ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx (we currently
    have a copy of this in data/pdb-20160711-entries.txt). Note that the file
    has TAB-separated fields.

    @param filename: A C{str} file name to read.
    @raises ValueError: If the input file is not in the expected format.
    @return: A C{dict} whose keys are C{str} PDB identifiers and whose
        values are C{dict}s with 'day', 'month', and 'year' keys with
        C{int} values.
    """

    result = {}

    headerLines = [
        'IDCODE, HEADER, ACCESSION DATE, COMPOUND, SOURCE, AUTHOR LIST, '
        'RESOLUTION, EXPERIMENT TYPE (IF NOT X-RAY)\n',
        '------- ------- --------------- --------- ------- ------------ '
        '----------- ---------------------------------------------------'
        '---------------------------------------------------------------'
        '---------------------------------------------------------------'
        '---------------------------------------------------------------'
        '----------------------------------------------------------\n'
    ]

    headerLineCount = len(headerLines)

    with open(filename) as fp:
        for lineNum, line in enumerate(fp, start=1):
            if lineNum <= headerLineCount:
                expected = headerLines[lineNum - 1]
                if expected is not None:
                    if line != expected:
                        raise ValueError(
                            'Line %d of %r was expected to be %r but was %r.'
                            % (lineNum, filename, expected, line))
            else:
                (pdbId, header, accessionDate, compound, source, authorList,
                 resolution, experimentType) = line.split('\t')

                if pdbId in result:
                    raise ValueError(
                        'Repeated PDB id %r found on line %d of input file %r.'
                        % (pdbId, lineNum, filename))

                month, day, year = map(int, accessionDate.split('/'))

                year += 1900 if year > 70 else 2000

                result[pdbId] = {
                    'day': day,
                    'month': month,
                    'year': year,
                }

    return result
