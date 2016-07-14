from __future__ import print_function

import sys
from json import dumps


def elmToJSON(classesFile, viralInstancesFile, outFile=sys.stdout):
    """
    Given a file with eukariotic linear motif (ELM) classes, and a file of
    viral ELM instances, this function will write out a line of JSON containing
    a regex 'pattern' and an ELM instance 'identifier'.

    Check data/README.md for a description of how to obtain the relevant files.

    @param classesFile: A C{str} filename of a file containing all ELM classes.
    @param viralInstancesFile: A C{str} filename of a file containing ELM
        instances specific to viruses.
    @param outFile: A file pointer where the output will be written to.
    """
    viralIdentifiers = set()

    lineCount = 0
    with open(viralInstancesFile) as vFp:
        for line in vFp:
            lineCount += 1
            if lineCount > 6:
                viralIdentifiers.add(line.split('\t')[2][1:-1])

    lineCount = 0
    with open(classesFile) as cFp:
        for line in cFp:
            lineCount += 1
            if lineCount > 6:
                splittedLine = line.split('\t')
                identifier = splittedLine[1][1:-1]
                pattern = splittedLine[4][1:-1]
                if identifier in viralIdentifiers:
                    print(dumps(
                        {
                            'identifier': identifier,
                            'pattern': pattern,
                        }, separators=(',', ':')), file=outFile)
