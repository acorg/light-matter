#!/usr/bin/env python

from operator import attrgetter

from light.landmarks import ALL_LANDMARK_FINDER_CLASSES
from light.trig import ALL_TRIG_FINDER_CLASSES

key = attrgetter('SYMBOL')

print 'Landmark finders:'
for finder in sorted(ALL_LANDMARK_FINDER_CLASSES, key=key):
    print '  %s = %s' % (finder.SYMBOL, finder.NAME)

print 'Trig point finders:'
for finder in sorted(ALL_TRIG_FINDER_CLASSES, key=key):
    print '  %s = %s' % (finder.SYMBOL, finder.NAME)
