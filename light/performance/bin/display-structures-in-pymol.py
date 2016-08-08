#!/usr/bin/env python

"""
A script which displays a set of given structures and their features in PyMOL.
"""

from __future__ import print_function

import argparse
from itertools import chain as ichain, repeat
from six.moves import zip
from Bio.PDB.PDBParser import PDBParser
from operator import attrgetter
import pymol
from pymol import cmd, stored

from dark.aa import find
from dark.reads import Reads, AAReadWithX

from light.backend import Backend
from light.database import Database
from light.landmarks import ALL_LANDMARK_CLASSES_INCLUDING_DEV
from light.parameters import DatabaseParameters, FindParameters
from light.trig import ALL_TRIG_CLASSES_INCLUDING_DEV
from light.performance.pymolgraphics import makeLegend

pymol.finish_launching()


# For a list of available colors and names, see:
# http://www.pymolwiki.org/index.php/Color_Values
COLORS = ['aquamarine', 'brightorange', 'darksalmon', 'deepsalmon',
          'deepolive', 'deepteal', 'firebrick', 'forest', 'greencyan',
          'lightblue', 'limon', 'marine', 'orange', 'palegreen', 'oxygen',
          'slate', 'tv_blue', 'tv_green', 'tv_red', 'violetpurple',
          'warmpink', 'violet', 'antimony', 'cesium', 'dubnium', 'gold',
          'iridium', 'phosphorus', 'magenta', 'yellow', 'zirconium',
          'actinium']

ALL_FEATURES = [
    (feature.SYMBOL, feature.NAME) for feature in
    sorted(ALL_LANDMARK_CLASSES_INCLUDING_DEV + ALL_TRIG_CLASSES_INCLUDING_DEV,
           key=attrgetter('NAME'))]

FEATURE_COLORS = dict(zip([feature[0] for feature in ALL_FEATURES], COLORS))

CHAIN_COLORS = ['grey', 'scandium', 'tin', 'deuterium', 'white']


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=('Displays a set of given structures and their features '
                     'in PyMOL.'))

    parser.add_argument(
        '--structureName', action='append', dest='structureNames',
        required=True,
        help=('The name of a structure to display. Can be specified multiple '
              'times. You can pass the same structure file multiple times '
              '(using --structureFile) but you will need to use a different '
              'structure name for each, otherwise the structures will be '
              'drawn on top of each other by PyMol. The structureName must be '
              'in the form of XXXX:Y, where X is the PDB identifier and Y is '
              'the chain that should be compared.'))

    parser.add_argument(
        '--structureFile', action='append', dest='structureFiles',
        required=True,
        help=('The structure file (in PDB format) to display. Must be '
              'specified the same number of times as --structureName '
              'is given.'))

    parser.add_argument(
        '--params', action='append',
        help=('Specify a full set of database and finder parameters (in '
              'quotes). This argument may be repeated. If less parameter '
              'sets are given than structure names, the final parameter '
              'set will be used repeatedly as many times as needed. If '
              'you need to put a space in an argument, use a "+". E.g., '
              '--params "--landmark PDB+AlphaHelix --landmark AC+AlphaHelix"'))

    parser.add_argument(
        '--colorBestBin', default=False, action='store_true',
        help=('If given, features in the best bin for the best matching '
              'chains of two structures will be colored red. Can only '
              'be used when exactly two structures are given.'))

    parser.add_argument(
        '--printParams', default=False, action='store_true',
        help='If given, print the values of all parameters used.')

    parser.add_argument(
        '--printColors', default=False, action='store_true',
        help='If given, print the colors for each finder.')

    parser.add_argument(
        '--showLegend', default=False, action='store_true',
        help=('If given, display a legend of which color corresponds to which '
              'feature.'))

    args = parser.parse_args()

    if args.colorBestBin and len(args.structureNames) != 2:
        raise ValueError(
            '--colorBestBin can only be used if 2 structures are given.')

    if len(args.structureFiles) != len(args.structureNames):
        raise ValueError('The number of structure files given (%d) must '
                         'match the number of structure names (%d).' %
                         (len(args.structureFiles), len(args.structureNames)))

    if args.params and len(args.params) > len(args.structureNames):
        raise ValueError('You cannot specify more parameter sets (%d) than '
                         'structure names (%d).' %
                         (len(args.params), len(args.structureNames)))

    if args.printColors:
        print('COLORS')
        for i, feature in enumerate(ALL_FEATURES):
            print(ALL_FEATURES[i][0], '\t', ALL_FEATURES[i][1], '\t',
                  COLORS[i])

    if args.params:
        givenParams = []
        for paramsStr in args.params:
            paramsParser = argparse.ArgumentParser()
            # Change '+' to ' ' to make it easy to specify finders that
            # have a space in their names.
            subArgs = [arg.replace('+', ' ') for arg in paramsStr.split()]
            FindParameters.addArgsToParser(paramsParser)
            DatabaseParameters.addArgsToParser(paramsParser)
            paramArgs = paramsParser.parse_args(subArgs)
            dbParams = DatabaseParameters.fromArgs(paramArgs)
            findParams = FindParameters.fromArgs(paramArgs)
            givenParams.append((dbParams, findParams))

        if len(givenParams) < len(args.structureNames):
            # We were given fewer parameter sets than structures.  Re-use
            # the final given parameters for structures that don't have a
            # correspondign set of given parameters.
            params = ichain(givenParams, repeat(givenParams[-1]))
        else:
            params = givenParams
    else:
        # The parameter settings for all structures will be the defaults.
        params = zip(repeat(DatabaseParameters()), repeat(FindParameters()))

    first = True

    # Loop over the structures that were provided and display them. We have
    # to use a version of zip that raises StopIteration when its shortest
    # iterator is consumed (params might have been set up in a simple way
    # above to iterate forever) so we use six.moves.zip which uses
    # itertools.izip when running under Python 2.
    for name, structureFile, (dbParams, findParams) in zip(
            args.structureNames, args.structureFiles, params):

        structureName, chainIdToCompare = name.split(':')

        if args.printParams:
            print('DATABASE PARAMETERS FOR STRUCTURE %r' % structureName)
            print(dbParams.print_(margin='  '))
            print('FIND PARAMETERS FOR STRUCTURE %r' % structureName)
            print(findParams.print_(margin='  '))

        # Set up the database.
        database = Database(dbParams)
        backend = Backend()
        backend.configure(database.dbParams)

        # Read the sequence out of the PDB file.
        s = PDBParser(PERMISSIVE=1).get_structure(structureName, structureFile)
        chains = Reads()
        for chain in s.get_chains():
            chainSequence = ''
            for aa in chain.get_residues():
                try:
                    aa1 = find(aa.resname).abbrev1
                except AttributeError:
                    aa1 = 'X'
                chainSequence += aa1
            chains.add(AAReadWithX(chain.id, chainSequence))

        # Load the structure into PyMOL.
        cmd.load(structureFile, structureName)

        # Set the display.
        cmd.show('cartoon')
        cmd.hide('lines')

        for chain in chains:
            print(chain.id)
            if chain.id == chainIdToCompare:
                chainToCompare = chain

        # Remember the name of the first structure. If we are going to
        # color the best bin, add its chains to its database. Align all
        # subsequent structures to the first structure.
        if first:
            first = False
            firstStructureName = structureName
            firstDatabase = database
            if args.colorBestBin:
                firstDatabase.addSubject(chainToCompare)
        else:
            cmd.align(firstStructureName, structureName)

        # Make sure the residues in each chain are zero-based.
        stored.first = None

        structureList = ['(model %s and (%s))' % (obj, structureName)
                         for obj in cmd.get_object_list(
                             '(' + structureName + ')')]
        structureChainList = ['(%s and chain %s)' % (structure, chain)
                              for structure in structureList for chain in
                              cmd.get_chains(structure)]

        for structureChain in structureChainList:
            cmd.iterate('first %s and polymer and n. CA' % structureChain,
                        'stored.first=resv')
            # Reassign the residue numbers.
            cmd.alter(structureChain,
                      'resi=str(int(resi)-%d)' % int(stored.first))
        cmd.rebuild()

        # Color the features and chains.
        for i, chain in enumerate(chains):
            # Color each chain.
            what = '%s & chain %s' % (structureName, chain.id)
            try:
                color = CHAIN_COLORS[i]
            except IndexError:
                color = 'white'
            cmd.color(color, what)

            # Color the features.
            scannedQuery = backend.scan(chain)
            for landmark in scannedQuery.landmarks:
                color = FEATURE_COLORS[landmark.symbol]
                start = landmark.offset
                end = landmark.offset + landmark.length
                what = 'resi %d-%d & %s & chain %s' % (start, end - 1,
                                                       structureName, chain.id)
                cmd.color(color, what)

            for trigPoint in scannedQuery.trigPoints:
                color = FEATURE_COLORS[trigPoint.symbol]
                start = trigPoint.offset
                end = trigPoint.offset + trigPoint.length
                what = 'resi %d-%d & %s & chain %s' % (start, end - 1,
                                                       structureName, chain.id)
                cmd.color(color, what)

        if args.colorBestBin and not first:
            result = firstDatabase.find(chainToCompare, findParams,
                                        storeFullAnalysis=True)
            analysis = result.analysis

            try:
                # Make sure there's no error if there's no significant bin
                significantBins = analysis['0']['significantBins']
            except KeyError:
                bin_ = {}
            else:
                bin_ = significantBins[0]['bin']
                subjectChain = \
                    firstDatabase.getSubjectByIndex('0').read

            for match in bin_:
                subjectLmStart = match['subjectLandmark'].offset
                subjectLmEnd = subjectLmStart + (
                    match['subjectLandmark'].length)
                subjectLandmark = 'resi %d-%d & %s & chain %s' % (
                    subjectLmStart, subjectLmEnd - 1, firstStructureName,
                    subjectChain.id)
                cmd.color('br9', subjectLandmark)

                subjectTpStart = match['subjectTrigPoint'].offset
                subjectTpEnd = subjectTpStart + (
                    match['subjectTrigPoint'].length)
                subjectTrigPoint = 'resi %d-%d & %s & chain %s' % (
                    subjectTpStart, subjectTpEnd - 1, firstStructureName,
                    subjectChain.id)
                cmd.color('br9', subjectTrigPoint)

                queryLmStart = match['queryLandmark'].offset
                queryLmEnd = queryLmStart + match['queryLandmark'].length
                queryLandmark = 'resi %d-%d & %s & chain %s' % (
                    queryLmStart, queryLmEnd - 1, structureName, chain.id)
                cmd.color('br9', queryLandmark)

                queryTpStart = match['queryTrigPoint'].offset
                queryTpEnd = queryTpStart + match['queryTrigPoint'].length
                queryTrigPoint = 'resi %d-%d & %s & chain %s' % (
                    queryTpStart, queryTpEnd - 1, structureName, chain.id)
                cmd.color('br9', queryTrigPoint)

    # Display structures in a grid.
    cmd.set('grid_mode')

    # Display the legend.
    if args.showLegend:
        cgo = makeLegend(ALL_FEATURES)
        cmd.load_cgo(cgo, 'legend')
        cmd.set('grid_slot', 0, 'legend')

    # Display the sequences.
    cmd.set('seq_view')
