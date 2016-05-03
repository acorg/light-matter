import sys

from dark.fasta import FastaReads
from dark.reads import AAReadWithX


def makeSubsets(length, helix):
    """
    Take a string and make subsets of a specified length.

    @param length: A C{int} of length of the subsets.
    @param helix: A {dark.AAReadsWithX} instance of sequence that should be
        subsetted.
    """
    assert(length < len(helix))
    numberOfHelices = len(helix) - length
    for i in range(numberOfHelices + 1):
        newHelixId = '%s[%d:%d]' % (helix.id, i, length + i)
        yield AAReadWithX(newHelixId, helix.sequence[i:length + i])


def makeAllSubsets(helixFile, outFile):
    """
    Take a fasta file of helices, and write out a fasta file with all possible
    helix subsets.

    @param helixFile: A C{str} filename of helices to be subsetted.
    @param outFile: A C{str} filename that the subsetted helices should be
        written to.
    """
    allHelices = [helix for helix in FastaReads(helixFile,
                                                readClass=AAReadWithX,
                                                checkAlphabet=0)]

    allLengths = [len(helix) for helix in allHelices]
    minLength = min(allLengths)
    maxLength = max(allLengths)

    with open(outFile, 'w') as ofp:
        for length in range(minLength, maxLength + 1):
            for helix in allHelices:
                if length == len(helix):
                    newHelix = AAReadWithX(helix.id, helix.sequence)
                    print(newHelix.toString(format_='fasta'), file=ofp)
                elif length < len(helix):
                    for newHelix in makeSubsets(length, helix):
                        print(newHelix.toString(format_='fasta'), file=ofp)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        from os.path import basename
        print('Usage: %s helixFile, outFile' % basename(sys.argv[0]),
              file=sys.stderr)
    else:
        helixFile, outFile = sys.argv[1:]
        makeAllSubsets(helixFile, outFile)
