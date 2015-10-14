#!/usr/bin/env python

import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.database import DatabaseComponent
from light.database import DatabaseSpecifier


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Start a WAMP-based distributed light-matter database.')

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    # We're always using WAMP for distributed databases.
    args.wampServer = True

    database = databaseSpecifier.getDatabaseFromArgs(args)
    runner = ApplicationRunner(args.wampUrl, args.realm,
                               extra=dict(database=database))
    runner.run(DatabaseComponent)