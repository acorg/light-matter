from dark.reads import AARead
from dark.fasta import FastaReads

from light.database import DatabaseSpecifier
from light.parameters import FindParameters


def affinityMatrix(sequences, findParams=None, **kwargs):
    """
    Produce an affinity matrix containing scores for a set of reads matched
    against the subjects in a database.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param findParams: A C{light.parameters.FindParameters} instance.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    @return: A two-dimensional array of match scores. The first dimension is
        the read number, the second is the database subject index.
    """
    if isinstance(sequences, str):
        reads = FastaReads(sequences, readClass=AARead, upperCase=True)
    else:
        reads = sequences

    findParams = findParams or FindParameters()

    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    subjectCount = db.subjectCount()
    affinity = []

    for readIndex, read in enumerate(reads):
        row = []
        append = row.append
        analysis = db.find(read, findParams).analysis
        for subjectIndex in map(str, range(subjectCount)):
            # Be careful how we access the analysis. It is a defaultdict,
            # so its keys are created upon access. I.e., use 'in' to test
            # for membership not try/except, because analysis[subjectIndex]
            # will never raise a KeyError.
            if subjectIndex in analysis:
                score = analysis[subjectIndex]['bestScore']
            else:
                score = 0.0
            append(score)

        affinity.append(row)

    return affinity
