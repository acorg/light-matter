#!/usr/bin/env python

import argparse
from operator import attrgetter

from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES

parser = argparse.ArgumentParser(
    description='List all landmark and trig point finders.')

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='If True, also print documentation for each finder class.')

args = parser.parse_args()
key = attrgetter('NAME')

maxLen = max(len(cls.NAME) for cls in
             set(ALL_LANDMARK_CLASSES) | set(ALL_TRIG_CLASSES))

print('%d landmark finders:' % len(ALL_LANDMARK_CLASSES))
for finder in sorted(ALL_LANDMARK_CLASSES, key=key):
    print('  %-*s (%s)' % (maxLen, finder.NAME, finder.SYMBOL))
    if args.verbose:
        print(finder.__doc__)

print('%d trig point finders:' % len(ALL_TRIG_CLASSES))
for finder in sorted(ALL_TRIG_CLASSES, key=key):
    print('  %-*s (%s)' % (maxLen, finder.NAME, finder.SYMBOL))
    if args.verbose:
        print(finder.__doc__)
