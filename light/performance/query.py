from collections import defaultdict

from dark.fasta import FastaReads


def queryDatabase(subjects, queries, database, findParams=None):
    """
    Add subjects to a database, query it, return results.

    @param subjects: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be turned into a database.
    @param queries: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be looked up in the database.
    @param database: A C{light.database.Database} instance.
    @param findParams: An instance of C{light.parameters.FindParameters} or
        C{None} to use default find parameters.
    @return: A C{dict} whose keys are query ids and whose values are C{dict}s
        that map subject ids to scores. I.e., for each read we provide a
        C{dict} showing what subjects it matched, and with what score.
    """
    if isinstance(queries, str):
        queries = FastaReads(queries, upperCase=True)

    if isinstance(subjects, str):
        subjects = FastaReads(subjects, upperCase=True)

    list(map(database.addSubject, subjects))

    resultDict = defaultdict(dict)

    for query in queries:
        result = database.find(query, findParams)
        for subjectIndex in result.significantSubjects():
            subject = database.getSubjectByIndex(subjectIndex)
            score = result.analysis[subjectIndex]['bestBinScore']
            resultDict[query.id][subject.id] = score

    return resultDict
