from collections import defaultdict

from dark.reads import Reads
from dark.fasta import FastaReads
from light.database import Database


def queryDatabase(databaseFile, queryFile, landmarkFinders, trigFinders,
                  limitPerLandmark=None, maxDistance=None, minDistance=None,
                  aboveMeanThreshold=None, bucketFactor=1):
    """
    Build a database, query it, return results.

    @param databaseFile: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be turned into a database.
    @param queryFile: Either an instance of C{dark.reads.Reads} or a C{str}
        filename with sequences that should be looked up in the database.
    @param landmarkFinders: a C{list} of landmarkFinders.
    @param landmarkFinders: a C{list} of landmarkFinders.
    @param limitPerLandmark: A limit on the number of pairs to yield per
        landmark per read.
    @param maxDistance: The maximum distance permitted between yielded pairs.
    @param minDistance: The minimum distance permitted between yielded pairs.
    @param bucketFactor: A C{int} factor by which the distance between
        landmark and trig point is divided, to influence sensitivity.
    @return: A C{dict} whose keys are query ids and whose values are C{dict}s
        that map subject ids to scores. I.e., for each read we provide a
        C{dict} showing what subjects it matched, and with what score.
    """
    database = Database(landmarkFinders, trigFinders, limitPerLandmark,
                        maxDistance, minDistance, bucketFactor)
    if isinstance(queryFile, Reads):
        queries = queryFile
    else:
        queries = FastaReads(queryFile)

    if isinstance(databaseFile, Reads):
        subjects = databaseFile
    else:
        subjects = FastaReads(databaseFile)
    map(database.addSubject, subjects)

    resultDict = defaultdict(dict)

    for query in queries:
        result = database.find(query, aboveMeanThreshold=aboveMeanThreshold)
        for subjectIndex in result.significant():
            subject = database.subjectInfo[subjectIndex][0]
            score = result.analysis[subjectIndex]['score']
            resultDict[query.id][subject] = score

    return resultDict
