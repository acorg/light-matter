#!/usr/bin/env python

"""
A script which displays a set of given structures and their features in PyMOL.
"""

from __future__ import print_function

import argparse
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
        help='The names of the structures that should be displayed')

    parser.add_argument(
        '--structureFile', action='append', dest='structureFiles',
        help=('The files of the structures that should be displayed in PDB '
              'file format'))

    parser.add_argument(
        '--colorBestBin', default=False, action='store_true',
        help=('If given, features in the best bin will be colored red. If '
              'there are multiple chains, then for each chain the best bin of '
              'the best match with any other chain will be colored.'))

    parser.add_argument(
        '--printParams', default=False, action='store_true',
        help='If given, print the values of all parameters used.')

    parser.add_argument(
        '--printColors', default=False, action='store_true',
        help='If given, print the colors for each finder.')

    parser.add_argument(
        '--showLegend', default=False, action='store_true',
        help='If given, display a legend of which color corresponds to which '
        'feature.')

    FindParameters.addArgsToParser(parser)
    DatabaseParameters.addArgsToParser(parser)

    args = parser.parse_args()

    dbParams = DatabaseParameters.fromArgs(args)
    findParams = FindParameters.fromArgs(args)

    if args.printParams:
        print('DATABASE PARAMETERS')
        print(dbParams.print_(margin='  '))
        print('FIND PARAMETERS')
        print(findParams.print_(margin='  '))

    if args.printColors:
        print('COLORS')
        for i, feature in enumerate(ALL_FEATURES):
            print(ALL_FEATURES[i][0], '\t', ALL_FEATURES[i][1], '\t',
                  COLORS[i])

    notFirst = False

    # Set up the database.
    database = Database(dbParams)
    backend = Backend()
    backend.configure(database.dbParams)

    # Loop over the structures that were provided and display them.
    for i, structureName in enumerate(args.structureNames):
        # Read the sequence out of the PDB file.
        s = PDBParser(PERMISSIVE=1).get_structure(structureName,
                                                  args.structureFiles[i])
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
        cmd.load(args.structureFiles[i])

        # Set the display.
        cmd.show('cartoon')
        cmd.hide('lines')

        # Align all structures to the first structure that was loaded.
        if not notFirst:
            notFirst = structureName
            if args.colorBestBin:
                for chain in chains:
                    database.addSubject(chain)
        else:
            cmd.align(notFirst, structureName)

        # Make sure the residues in each chain are zero-based.
        stored.first = None

        structureList = ['(model %s and (%s))' % (p, structureName)
                         for p in cmd.get_object_list(
            '(' + structureName + ')')]
        structureChainList = ['(%s and chain %s)' % (structure, chain)
                              for structure in structureList for chain in
                              cmd.get_chains(structure)]

        for structureChain in structureChainList:
            cmd.iterate('first %s and polymer and n. CA' % structureChain,
                        'stored.first=resv')
            # reassign the residue numbers.
            cmd.alter('%s' % structureChain,
                      'resi=str(int(resi)-%s)' % str(int(stored.first)))
        cmd.rebuild()

        # Color the features and chains.
        for i, chain in enumerate(chains):
            # Color each chain in a shade of grey.
            what = '%s & chain %s' % (structureName, chain.id)
            try:
                color = CHAIN_COLORS[i]
            except IndexError:
                color = 'white'
            cmd.color(color, what)

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
                what = 'resi %d-%d & chain %s & %s' % (start, end - 1,
                                                       chain.id, structureName)
                cmd.color(color, what)

    if args.colorBestBin:
        if notFirst != structureName:
            for chain in chains:
                result = database.find(chain, findParams,
                                       storeFullAnalysis=True)
                analysis = result.analysis
                # Get the best match
                bestScore = 0
                bestSbjctInd = 0
                for subjectIndex in analysis:
                    if analysis[subjectIndex]['overallScore'] > bestScore:
                        bestScore = analysis[subjectIndex]['overallScore']
                        bestSbjctInd = subjectIndex
                try:
                    # Make sure there's no error if there's no significant bin
                    significantBins = analysis[bestSbjctInd]['significantBins']
                except KeyError:
                    bin_ = {}
                else:
                    bin_ = significantBins[0]['bin']
                    subjectChain = \
                        database.getSubjectByIndex(bestSbjctInd).read

                for match in bin_:
                    subjectLmStart = match['subjectLandmark'].offset
                    subjectLmEnd = subjectLmStart + (
                        match['subjectLandmark'].length)
                    subjectLandmark = 'resi %d-%d & %s & chain %s' % (
                        subjectLmStart, subjectLmEnd - 1, notFirst,
                        subjectChain.id)
                    cmd.color('br9', subjectLandmark)

                    subjectTpStart = match['subjectTrigPoint'].offset
                    subjectTpEnd = subjectTpStart + (
                        match['subjectTrigPoint'].length)
                    subjectTrigPoint = 'resi %d-%d & %s & chain %s' % (
                        subjectTpStart, subjectTpEnd - 1, notFirst,
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

    # Display the structures next to each other.
    cmd.set('grid_mode')

    # Display the legend.
    if args.showLegend:
        cgo = makeLegend(ALL_FEATURES)
        cmd.load_cgo(cgo, 'legend')
        cmd.set('grid_slot', 0, 'legend')

    # Display the sequences.
    cmd.set('seq_view')
