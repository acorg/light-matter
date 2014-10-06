import MySQLdb
from os import environ


def getDatabaseConnection():
    db = MySQLdb.connect(host='localhost',
                         user=environ.get('DBI_USER', environ['USER']),
                         passwd=environ['DBI_PASSWORD'],
                         db='lightmatter')
    return db
