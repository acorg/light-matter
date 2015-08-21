#!/usr/bin/env python

import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.shutdown import ShutdownComponent
from light.database import DatabaseSpecifier


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Shut down (and, optionally, save) a WAMP-based '
                     'distributed light-matter database.'))

    parser.add_argument(
        '--noSave', action='store_true', default=False,
        help='If True, shut down the database without saving.')

    databaseSpecifier = DatabaseSpecifier(allowFromFile=False,
                                          allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    args.wamp = True  # We're always using WAMP for distributed databases.

    runner = ApplicationRunner(args.wampUrl, args.realm,
                               extra=dict(noSave=args.noSave,
                                          filePrefix=args.filePrefix))
    runner.run(ShutdownComponent)
