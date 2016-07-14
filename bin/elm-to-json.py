#!/usr/bin/env python

from __future__ import print_function

import sys

from light.elm import elmToJSON

if len(sys.argv) != 3:
    print('Usage: %s elm-classes elm-viral-instances' % sys.argv[0],
          file=sys.stderr)
    sys.exit(1)
else:
    classesFile = sys.argv[1]
    viralInstancesFile = sys.argv[2]

    elmToJSON(classesFile, viralInstancesFile)
