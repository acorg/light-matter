#!/usr/bin/env python

"""
Read the PDB ss.txt file on stdin and do the following:

  1) Change sequence id lines like ">101M:A:sequence" to ">pdb_101m_a"
  2) Change structure id lines like ">101M:A:secstr" to ">pdb_101m_a:structure"
  3) Replace spaces (present in structure) with hyphens
  4) Check each sequence id is followed by a matching structure id

Write the modified output to stdout.

The input will typically be the ss.txt secondary structure file available at
http://www.rcsb.org/pdb/static.do?p=download/http/index.html
"""

import sys

from light.performance.utils import pdbNameToPythonName

for line in sys.stdin:
    line = line[:-1]
    if line[0] == '>':
        if line.endswith(':sequence'):
            sequenceName = line[1:-9]
            print('>%s' % pdbNameToPythonName(sequenceName))
        elif line.endswith(':secstr'):
            structureName = line[1:-7]
            # Make sure the name in the structure part matches the name we
            # saw in the sequence section.
            assert structureName == sequenceName, (
                "Sequence name %r doesn't match the following structure "
                "'name %r." % (sequenceName, structureName))
            print('>%s:structure' % pdbNameToPythonName(structureName))
        else:
            raise RuntimeError('Unrecognized sequence id line %r' % line)
    else:
        print(line.replace(' ', '-'))
