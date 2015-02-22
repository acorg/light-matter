from collections import defaultdict

from dark.fasta import FastaReads


def queryDatabase(subjects, queries, database, significanceFraction=None):
    """
    Add subjects to a database, query it, return results.

    @param subjects: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be turned into a database.
    @param queries: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be looked up in the database.
    @param database: A C{light.database.Database} instance.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @return: A C{dict} whose keys are query ids and whose values are C{dict}s
        that map subject ids to scores. I.e., for each read we provide a
        C{dict} showing what subjects it matched, and with what score.
    """
    if isinstance(queries, basestring):
        queries = FastaReads(queries)

    if isinstance(subjects, basestring):
        subjects = FastaReads(subjects)

    map(database.addSubject, subjects)

    resultDict = defaultdict(dict)

    for query in queries:
        result = database.find(query,
                               significanceFraction=significanceFraction)
        for subjectIndex in result.significant():
            subject = database.subjectInfo[subjectIndex][0]
            score = result.analysis[subjectIndex]['score']
            resultDict[query.id][subject] = score

    return resultDict
