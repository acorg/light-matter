#!/usr/bin/env python

"""
A script which displays a set of given structures and their features in PyMOL.
"""

from __future__ import print_function

import argparse
from itertools import chain as ichain, repeat
from six.moves import zip
from operator import attrgetter
from os.path import join, dirname
import pymol
from pymol import cmd, stored

from dark.reads import Reads, SSAAReadWithX
from dark.fasta_ss import SSFastaReads

import light
from light.backend import Backend
from light.database import Database
from light.landmarks import ALL_LANDMARK_CLASSES_INCLUDING_DEV
from light.parameters import DatabaseParameters, FindParameters
from light.trig import ALL_TRIG_CLASSES_INCLUDING_DEV
from light.performance.pymolgraphics import makeLegend

from pymol.cgo import VERTEX, END, LINES, COLOR, BEGIN

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
        '--showPairs', default=False, action='store_true',
        help=('If given, lines will be drawn between matching features in the '
              'best bin. Can only be used together with colorBestBin.'))

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

    if (len(args.structureNames) == 2 and args.colorBestBin is True and
            args.showPairs is True):
        setGridMode = False
    else:
        setGridMode = True

    # Loop over the structures that were provided and display them. We have
    # to use a version of zip that raises StopIteration when its shortest
    # iterator is consumed (params might have been set up in a simple way
    # above to iterate forever) so we use six.moves.zip which uses
    # itertools.izip when running under Python 2.
    for name, structureFile, (dbParams, findParams) in zip(
            args.structureNames, args.structureFiles, params):

        structureName, chainIdToCompare = name.split(':')
        chainIdToCompare = chainIdToCompare.lower()

        if args.printParams:
            print('DATABASE PARAMETERS FOR STRUCTURE %r' % structureName)
            print(dbParams.print_(margin='  '))
            print('FIND PARAMETERS FOR STRUCTURE %r' % structureName)
            print(findParams.print_(margin='  '))

        # Set up the database.
        database = Database(dbParams)
        backend = Backend()
        backend.configure(database.dbParams)

        # Read the sequence out of the PDB ss.txt file.
        chains = Reads()
        sequenceFile = join(dirname(light.__file__),
                            '..', 'data', 'pdb-20160303-ss.txt')
        for record in SSFastaReads(sequenceFile, checkAlphabet=0):
            if structureName in record.id:
                chainName = record.id.split('_')[2].lower()
                chains.add(SSAAReadWithX(chainName, record.sequence,
                                         record.structure))

        assert len(chains) > 0, ('%r does not contain any sequences with id %r'
                                 % (sequenceFile, structureName))
        if first:
            firstStructureName = structureName + '1'
            structureName += '1'

        # Load the structure into PyMOL.
        cmd.load(structureFile, structureName)

        # Set the display.
        cmd.show('cartoon')
        cmd.hide('lines')

        # Keep a record of the chain that should be compared.
        for chain in chains:
            if chain.id == chainIdToCompare:
                chainToCompare = chain
                break
        if not chainToCompare:
            raise ValueError('%r has no chain with name %r.' % (
                structureName, chainIdToCompare))

        # Remember the name of the first structure. If we are going to
        # color the best bin, add its chains to its database. Align all
        # subsequent structures to the first structure.
        if first:
            firstChainToCompare = chainToCompare
            firstDatabase = database
            if args.colorBestBin:
                firstDatabase.addSubject(chainToCompare)
        else:
            cmd.align('%s & chain %s' % (firstStructureName,
                                         firstChainToCompare.id.upper()),
                      '%s & chain %s' % (structureName,
                                         chainToCompare.id.upper()))

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

        # Alter the coordinates if necessary.
        if first:
            shift = 100
        else:
            shift = -100
        if not setGridMode:
            print('altering', structureName, shift)
            cmd.alter_state(1, structureName, 'x+=%d' % shift)

        cmd.rebuild()

        # Color the features and chains.
        for i, chain in enumerate(chains):
            if chain.id == chainToCompare.id:
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
                                                           structureName,
                                                           chain.id)
                    cmd.color(color, what)

                for trigPoint in scannedQuery.trigPoints:
                    color = FEATURE_COLORS[trigPoint.symbol]
                    start = trigPoint.offset
                    end = trigPoint.offset + trigPoint.length
                    what = 'resi %d-%d & %s & chain %s' % (start, end - 1,
                                                           structureName,
                                                           chain.id)
                    cmd.color(color, what)
            else:
                what = '%s & chain %s' % (structureName, chain.id.upper())
                cmd.color('black', what)

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
            print(len(bin_))
            cgo = [
                BEGIN, LINES,
                COLOR, 0.8, 0.8, 0.8]
            for match in bin_:
                subjectLmStart = match['subjectLandmark'].offset
                subjectLmEnd = subjectLmStart + (
                    match['subjectLandmark'].length)
                subjectLandmark = 'resi %d-%d & %s & chain %s' % (
                    subjectLmStart, subjectLmEnd - 1, firstStructureName,
                    subjectChain.id)
                cmd.color('br9', subjectLandmark)
                slmxyz = cmd.get_model(
                    'n. CA & resi %d, & %s & chain %s' % (
                        subjectLmStart, firstStructureName, subjectChain.id),
                    1).get_coord_list()

                subjectTpStart = match['subjectTrigPoint'].offset
                subjectTpEnd = subjectTpStart + (
                    match['subjectTrigPoint'].length)
                subjectTrigPoint = 'resi %d-%d & %s & chain %s' % (
                    subjectTpStart, subjectTpEnd - 1, firstStructureName,
                    subjectChain.id)
                cmd.color('br9', subjectTrigPoint)
                stpxyz = cmd.get_model(
                    'n. CA & resi %d, & %s & chain %s' % (
                        subjectTpStart, firstStructureName, subjectChain.id),
                    1).get_coord_list()

                queryLmStart = match['queryLandmark'].offset
                queryLmEnd = queryLmStart + match['queryLandmark'].length
                queryLandmark = 'resi %d-%d & %s & chain %s' % (
                    queryLmStart, queryLmEnd - 1, structureName,
                    chainToCompare.id)
                cmd.color('br9', queryLandmark)
                qlmxyz = cmd.get_model(
                    'n. CA & resi %d, & %s & chain %s' % (
                        queryLmStart, structureName,
                        chainToCompare.id), 1).get_coord_list()

                queryTpStart = match['queryTrigPoint'].offset
                queryTpEnd = queryTpStart + match['queryTrigPoint'].length
                queryTrigPoint = 'resi %d-%d & %s & chain %s' % (
                    queryTpStart, queryTpEnd - 1, structureName,
                    chainToCompare.id)
                cmd.color('br9', queryTrigPoint)
                qtpxyz = cmd.get_model(
                    'n. CA & resi %d, & %s & chain %s' % (
                        queryTpStart, structureName,
                        chainToCompare.id), 1).get_coord_list()

                try:
                    cgo.extend([VERTEX, slmxyz[0][0], slmxyz[0][1],
                                slmxyz[0][2],
                                VERTEX, qlmxyz[0][0], qlmxyz[0][1],
                                qlmxyz[0][2]])
                except IndexError:
                    continue

                try:
                    cgo.extend([VERTEX, stpxyz[0][0], stpxyz[0][1],
                                stpxyz[0][2],
                                VERTEX, qtpxyz[0][0], qtpxyz[0][1],
                                qtpxyz[0][2]])
                except IndexError:
                    continue

        if first:
            first = False

    if not setGridMode:
        cgo.append(END)
        cmd.load_cgo(cgo, 'pairs')

    # Display structures in a grid.
    if setGridMode:
        cmd.set('grid_mode')

    # Display the legend.
    if args.showLegend:
        cgo = makeLegend(ALL_FEATURES)
        cmd.load_cgo(cgo, 'legend')
        cmd.set('grid_slot', 0, 'legend')

    # Display the sequences.
    cmd.set('seq_view')

    print('Done')
