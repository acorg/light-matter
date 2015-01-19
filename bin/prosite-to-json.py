#!/usr/bin/env python

import sys

from light.prosite import prositeToJSON

if len(sys.argv) != 2:
    print >>sys.stderr, 'Usage: %s prosite.dat' % sys.argv[0]
    sys.exit(1)
else:
    filename = sys.argv[1]

    prositeToJSON(filename)
