#!/usr/bin/env python

import argparse

from dark.reads import AARead

from light.database import DatabaseSpecifier


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Scan a sequence for landmarks and trig points')

    parser.add_argument(
        'sequence', help='Amino acid sequence to scan')

    databaseSpecifier = DatabaseSpecifier(allowFromFile=False,
                                          allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    database = databaseSpecifier.getDatabaseFromArgs(args)

    scannedRead = database.scan(AARead('id', args.sequence))

    if scannedRead.landmarks:
        print('All %d landmarks:' % len(scannedRead.landmarks))
        for landmark in scannedRead.landmarks:
            print('  ', landmark)
    else:
        print('No landmarks.')

    print()
    if scannedRead.trigPoints:
        print('All %d trig points:' % len(scannedRead.trigPoints))
        for trigPoint in scannedRead.trigPoints:
            print('  ', trigPoint)
    else:
        print('No trig points.')

    print()
    hashes = database.getHashes(scannedRead)
    if hashes:
        print('%s hash%s (according to database parameters):' % (
            len(hashes), '' if len(hashes) == 1 else 'es'), end='')
        for hash_, hashInfo in hashes.items():
            print()
            print('  ', hashInfo['landmark'])
            print('  ', hashInfo['trigPoint'])
            print('  ', hash_)
    else:
        print('No Hashes (according to database parameters).')
