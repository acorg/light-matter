import numpy as np

from dark.reads import AARead
from dark.fasta import FastaReads

from light.database import DatabaseSpecifier
from light.parameters import FindParameters


def affinityMatrix(queries, findParams=None, subjects=None, symmetric=True,
                   **kwargs):
    """
    Produce an affinity matrix containing scores for a set of reads matched
    against the subjects in a database.

    @param queries: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param findParams: A C{light.parameters.FindParameters} instance.
    @param subjects: Either 1) a C{str} filename of sequences to consider, or
        2) a C{light.reads.Reads} instance, or 3) C{None} to use the
        C{queries} as the subjects.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    @return: A two-dimensional array of match scores. The first dimension is
        the read number, the second is the database subject index.
    """
    if isinstance(queries, str):
        queries = list(FastaReads(queries, readClass=AARead, upperCase=True))

    if subjects is None:
        subjects = queries
    else:
        if isinstance(subjects, str):
            subjects = list(FastaReads(subjects, readClass=AARead,
                                       upperCase=True))

    findParams = findParams or FindParameters()

    db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)

    subjectIndices = []
    for subject in subjects:
        _, subjectIndex, _ = db.addSubject(subject)
        subjectIndices.append(subjectIndex)

    nQueries = len(queries)
    nSubjects = len(subjectIndices)

    affinity = np.zeros((nQueries, nSubjects))

    for i, query in enumerate(queries):
        analysis = db.find(query, findParams).analysis
        for j in range(nSubjects):
            if j < i and symmetric:
                score = affinity[j][i]
            else:
                # Be careful how we access the analysis. It is a defaultdict,
                # so its keys are created upon access. I.e., use 'in' to test
                # for membership not try/except, because analysis[subjectIndex]
                # will never raise a KeyError.
                if subjectIndices[j] in analysis:
                    score = analysis[subjectIndices[j]]['bestScore']
                else:
                    score = 0.0
            affinity[i][j] = score

    return affinity
