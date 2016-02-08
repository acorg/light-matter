#!/usr/bin/env python

from __future__ import print_function

import argparse

from light.landmarks import ALL_LANDMARK_CLASSES_EVEN_BAD_ONES
from light.trig import ALL_TRIG_CLASSES_EVEN_BAD_ONES


parser = argparse.ArgumentParser(
    description='List all landmark and trig point finders.')

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='If True, also print documentation for each finder class.')

args = parser.parse_args()

finders = ALL_LANDMARK_CLASSES_EVEN_BAD_ONES
print('%d landmark finders:' % len(finders))
for finder in finders:
    print('  ', finder.NAME)
    if args.verbose:
        print('  ', finder.__doc__)

finders = ALL_TRIG_CLASSES_EVEN_BAD_ONES
print('%d trig point finders:' % len(finders))
for finder in finders:
    print('  ', finder.NAME)
    if args.verbose:
        print('  ', finder.__doc__)
