#!/usr/bin/env python

import argparse

from dark.reads import AARead

from light.database import DatabaseSpecifier


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Examine an amino acid sequence for its landmarks and '
                     'trig points, and print them all, as well as the hashes '
                     'that would be used by the database. Note that ALL '
                     'landmarks and trig points are printed, but the '
                     'particular database parameters might mean that not all '
                     'of them are used in practice, depending on parameters '
                     'like maxDistance, limitPerLandmark, etc.'))

    parser.add_argument(
        'sequence', help='The amino acid sequence to examine')

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
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
        print('%d hash%s (according to database parameters):' % (
            len(hashes), '' if len(hashes) == 1 else 'es'), end='')
        for hash_, hashInfo in hashes.items():
            print()
            print('  ', hashInfo['landmark'])
            print('  ', hashInfo['trigPoint'])
            print('  ', hash_)
    else:
        print('No Hashes (according to database parameters).')
