#!/usr/bin/env python

from __future__ import print_function

import argparse

from light.landmarks import ALL_LANDMARK_CLASSES_EVEN_BAD_ONES
from light.trig import ALL_TRIG_CLASSES_EVEN_BAD_ONES
from light.parameters import DatabaseParameters


parser = argparse.ArgumentParser(
    description='List all landmark and trig point finders.')

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='If True, also print documentation for each finder class.')

args = parser.parse_args()
dbParams = DatabaseParameters(
    landmarks=ALL_LANDMARK_CLASSES_EVEN_BAD_ONES,
    trigPoints=ALL_TRIG_CLASSES_EVEN_BAD_ONES)

finders = dbParams.landmarkFinders
print('%d landmark finders:' % len(finders))
for finder in finders:
    print(finder.print_(margin='  '))
    if args.verbose:
        print('  ', finder.__doc__)

finders = dbParams.trigPointFinders
print('%d trig point finders:' % len(finders))
for finder in finders:
    print(finder.print_(margin='  '))
    if args.verbose:
        print('  ', finder.__doc__)
