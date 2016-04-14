import numpy as np

from dark.reads import AAReadWithX
from dark.fasta import FastaReads

from light.database import DatabaseSpecifier
from light.parameters import FindParameters


def affinityMatrix(queries, findParams=None, subjects=None, symmetric=True,
                   computeDiagonal=False, diagonalValue=1.0, progressFunc=None,
                   returnDict=False, **kwargs):
    """
    Produce an affinity matrix containing scores for a set of reads matched
    against a set of subjects.

    @param queries: Either A C{str} filename of sequences to consider or
        a C{light.reads.Reads} instance.
    @param findParams: A C{light.parameters.FindParameters} instance.
    @param subjects: Either 1) a C{str} filename of sequences to consider, or
        2) a C{light.reads.Reads} instance, or 3) C{None}, in which case
        the C{queries} will also be used as the subjects.
    @param symmetric: If C{True}, corresponding off-diagonal scores will be
        assumed to be equal and only computed once. I.e., this is a speedup
        when scores (affinities) are symmetric. This option gives the biggest
        speed up on a square matrix, but can also be used when the matrix is
        not square (e.g., when making a 2x4 matrix comparing {A, B} against
        {A, B, C, D}, the A->B distance can be used to set the B->A distance).
    @param computeDiagonal: If C{True}, values on the diagonal will be computed
        (i.e., obtained from find). Otherwise, all diagonal values will be set
        to C{diagonalValue}.
    @param diagonalValue: The result that diagonal values will all be set to if
        C{computeDiagonal} is False.
    @param progressFunc: If not C{None}, a function that takes two arguments.
        The function will be called before each query sequence is processed.
        The arguments will be the C{int} (zero-based) number of the query and
        the query (an AAReadWithX instance) itself.
    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    @param returnDict: If C{True}, return a C{dict} keyed by query id, whose
        values are C{dict}s keyed by subject id, whose values are C{float}
        scores. In other words, a 2-level deep C{dict} that allows the caller
        to look up a score via something like C{result[query.id][subject.id]}.
    @raise ValueError: If C{returnDict} is C{True} and there is a duplicated
        query or subject id.
    @return: If C{returnDict} is C{True}, a C{dict} as just described, else a
        two-dimensional array of match scores. The first dimension is the query
        index (in C{queries}, the second is the subject index (in C{subjects}).
    """
    if isinstance(queries, str):
        queries = list(FastaReads(queries, readClass=AAReadWithX,
                       upperCase=True))

    if subjects is None:
        subjects = queries
    else:
        if isinstance(subjects, str):
            subjects = list(FastaReads(subjects, readClass=AAReadWithX,
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
        if progressFunc:
            progressFunc(i, query)
        if symmetric:
            # We don't have to consider all subjects in the find, so pass a
            # restricted set of subject indices to restrict the search to.
            # The ones we omit have already been looked up.
            #
            # For clarity, there's a little code repetition here.
            if computeDiagonal:
                wantedIndices = set(subjectIndices[i:])
            else:
                wantedIndices = set(subjectIndices[i + 1:])
            result = db.find(query, findParams, subjectIndices=wantedIndices)
        else:
            result = db.find(query, findParams)

        analysis = result.analysis
        for j in range(nSubjects):
            if j < i and symmetric:
                score = affinity[j][i]
            elif j == i and not computeDiagonal:
                score = diagonalValue
            else:
                # Be careful how we access the analysis. It is a defaultdict,
                # so its keys are created upon access. I.e., use 'in' to test
                # for membership not try/except, because analysis[subjectIndex]
                # will never raise a KeyError.
                if subjectIndices[j] in analysis:
                    score = analysis[subjectIndices[j]]['bestBinScore']
                else:
                    score = 0.0
            affinity[i][j] = score

    if returnDict:
        result = {}
        for i, query in enumerate(queries):
            if query.id in result:
                raise ValueError('Query id %r appears more than once.' %
                                 query.id)
            result[query.id] = values = {}
            for j, subject in enumerate(subjects):
                if subject.id in values:
                    raise ValueError('Subject id %r appears more than once.' %
                                     subject.id)
                values[subject.id] = affinity[i][j]
        return result
    else:
        return affinity
