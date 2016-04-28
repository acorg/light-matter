#!/usr/bin/env python

"""
Print all database and find parameters available via the command line,
with their default values (according to argparse).
"""

import argparse

from light.parameters import DatabaseParameters, FindParameters

dbParams = DatabaseParameters()
parser = argparse.ArgumentParser(
    add_help=False,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
dbParams.addArgsToParser(parser)
print('Database parameters:')
# Splitting on 'optional arguments:\n' is a hack that works nicely as the
# first part of the help message is a summary of how to call the program
# and doesn't show descriptions or default values.
print(parser.format_help().split('optional arguments:\n')[1])

findParams = FindParameters()
parser = argparse.ArgumentParser(
    add_help=False,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
findParams.addArgsToParser(parser)
print('Find parameters:')
print(parser.format_help().split('optional arguments:\n')[1])
