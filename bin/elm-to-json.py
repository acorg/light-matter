#!/usr/bin/env python

from __future__ import print_function

import sys

from light.elm import elmToJSON

if len(sys.argv) == 3:
    classesFile = sys.argv[1]
    viralInstancesFile = sys.argv[2]
    elmToJSON(classesFile, viralInstancesFile=viralInstancesFile)
elif len(sys.argv) == 2:
    classesFile = sys.argv[1]
    elmToJSON(classesFile, viralInstancesFile=None)
else:
    print('Usage: %s elm-classes elm-viral-instances' % sys.argv[0],
          file=sys.stderr)
    sys.exit(1)
