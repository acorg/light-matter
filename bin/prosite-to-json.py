#!/usr/bin/env python

from __future__ import print_function

import sys

from light.prosite import prositeToJSON

if len(sys.argv) != 2:
    print('Usage: %s prosite.dat' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    filename = sys.argv[1]

    prositeToJSON(filename)
