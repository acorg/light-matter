#!/usr/bin/env python

import sys
import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.backend import BackendComponent
from light.database import DatabaseSpecifier

if sys.version_info < (3, 3):
    raise Exception('The light matter autobahn code needs Python 3.3 or '
                    'later.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Start a WAMP-based distributed light-matter database '
                     'backend.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()

    runner = ApplicationRunner(args.wampUrl, args.realm,
                               extra=dict(args=args))
    runner.run(BackendComponent)
