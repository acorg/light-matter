from dark.reads import AARead
from dark.fasta import FastaReads

from light.database import DatabaseSpecifier


def affinityMatrix(sequences, **kwargs):
    """
    Produce an affinity matrix containing scores for a set of reads matched
    against the subjects in a database.

    @param sequences: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    @return: A two-dimensional array of match scores. The first dimension is
        the read number, the second is the database subject index.
    """
    if isinstance(sequences, basestring):
        reads = FastaReads(sequences, readClass=AARead)
    else:
        reads = sequences

    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
    subjectCount = db.subjectCount
    affinity = []

    for readIndex, read in enumerate(reads):
        row = []
        append = row.append
        analysis = db.find(read).analysis
        for subjectIndex in xrange(subjectCount):
            try:
                score = analysis[subjectIndex]['score']
            except KeyError:
                score = 0.0
            append(score)

        affinity.append(row)

    return affinity
