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
        # The logic here may seem backward. We want the default to be that
        # databases are saved when shut down. The user shouldn't have to
        # specify anything on the command line to get the default behavior.
        # So in order to not save, the user will specify --noSave. But in
        # code, it's simpler to work with positive variable names, so we'd
        # like to write 'if save' rather than needing double negatives like
        # 'if not noSave' or 'noSave=False' (in method definitions).
        #
        # To get both these things at once, we present the user with an
        # option called '--noSave' but we store (dest) the command line
        # value into the 'save' attribute (whose value is the opposite of
        # noSave). For this reason, the default of the variable is True and
        # the argparse action is to save a False value. That makes for a
        # little awkwardness here (and this long comment) but everywhere
        # else in our code we just have simple 'save' variables that
        # default to True and we avoid all the double negatives.
        '--noSave', action='store_false', default=True, dest='save',
        help='If True, shut down the database WITHOUT saving.')

    databaseSpecifier = DatabaseSpecifier(allowCreation=False,
                                          allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()

    # We know we're a WAMP client, so don't needlessly make the user
    # specify --wampClient.
    args.wampClient = True

    runner = ApplicationRunner(
        args.wampUrl, args.realm,
        extra=dict(save=args.save, filePrefix=args.filePrefix))
    runner.run(ShutdownComponent)
