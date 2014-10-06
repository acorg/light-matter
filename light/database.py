from light.statistics.alpha_helix import AlphaHelix
from light import mysql
from dark.reads import Reads


def addToDatabase(fastaFile, statistic=AlphaHelix, db=None):
    """
    A script that adds calculated statistics to database.

    @param fastaFile: a C{str} title of a fasta file containing the sequences
        for which the statistic should be calculated.
    @param statistic: the statistic for which row should be added to database.
    @param db: A database connection, with C{cursor} and C{close} methods.
    """
    # open database connection, create table.
    if db is None:
        openedDb = True
        db = mysql.getDatabaseConnection()
    else:
        openedDb = False
    cursor = db.cursor()
    createTable = 'CREATE TABLE %s (virus VARCHAR(300), '
    'AlphaHelix VARCHAR(300));' % str(statistic)
    cursor.execute(createTable)

    # calculate statistic, add to table.
    stat = statistic()

    viruses = Reads(fastaFile)

    for virus in viruses:
        dist = stat._evaluate(virus, distances=True)
        query = 'INSERT INTO %s (virus, AlphaHelix) '
        'VALUES (%s, %s)' % (str(statistic), str(virus.id), str(dist))
        cursor.execute(query)

    if openedDb:
        db.close()


def compare(readStat, statistic=AlphaHelix, db=None):
    """
    A function that compares the values in the database with the values
    calculated for a read.

    @param readStat: a statistic calculated for a read.
    @param statistic: the statistic for which row should be added to database.
    """
    # open database connection.
    if db is None:
        openedDb = True
        db = mysql.getDatabaseConnection()
    else:
        openedDb = False
    cursor = db.cursor()

    # look in database for helix motif.
    query = 'SELECT virus FROM %s WHERE %s = %s' % (str(statistic),
                                                    str(statistic),
                                                    str(readStat))
    cursor.execute(query)
    presentVirus = cursor.fetchone()

    if openedDb:
        db.close()

    return presentVirus
