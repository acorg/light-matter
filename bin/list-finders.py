#!/usr/bin/env python

import argparse
from operator import attrgetter

from light.landmarks import ALL_LANDMARK_FINDER_CLASSES
from light.trig import ALL_TRIG_FINDER_CLASSES

parser = argparse.ArgumentParser(
    description='List all landmark and trig point finders.')

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='If True, also print documentation for each finder class.')

args = parser.parse_args()
key = attrgetter('NAME')

print 'Landmark finders:'
for finder in sorted(ALL_LANDMARK_FINDER_CLASSES, key=key):
    print '  %s (symbol %s)' % (finder.NAME, finder.SYMBOL)
    if args.verbose:
        print finder.__doc__

print 'Trig point finders:'
for finder in sorted(ALL_TRIG_FINDER_CLASSES, key=key):
    print '  %s (symbol %s)' % (finder.NAME, finder.SYMBOL)
    if args.verbose:
        print finder.__doc__
