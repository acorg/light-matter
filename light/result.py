class ScannedReadDatabaseResult(object):
    """
    A class that holds the results from a database lookup.

    @param subjectId: a C{str} name of the subject that was matched.
    @param queryId: a C{str} name of the query that was matched.
    @param offset: the difference in the offset of the feature in the subject
        minus the difference in the offset of the query
    """
    def __init__(self, subjectId, queryId, offset):
        self.subjectId = subjectId
        self.queryId = queryId
        self.offset = offset
