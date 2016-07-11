#!/usr/bin/env python

"""
A script which displays a set of given structures and their features in PyMOL.
"""

from __future__ import print_function

from Bio.PDB.PDBParser import PDBParser
import argparse
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


COLORS = ['aquamarine', 'brightorange', 'darksalmon', 'deepolive',
          'deepsalmon', 'deepteal', 'firebrick', 'forest', 'greencyan',
          'lightblue', 'limon', 'marine', 'orange', 'palegreen', 'oxygen',
          'slate', 'tv_blue', 'tv_green', 'tv_red', 'violetpurple',
          'warmpink', 'violet', 'antimony', 'cesium', 'dubnium', 'gold',
          'iridium', 'phosphorus']

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
        '--printParams', default=False, action='store_true',
        help='If given, print the values of all parameters used.')

    parser.add_argument(
        '--printColors', default=False, action='store_true',
        help='If given, print the colors for each finder.')

    parser.add_argument(
        '--makeLegend', default=False, action='store_true',
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

    if args.makeLegend:
        cgo = makeLegend()
        cmd.load_cgo(cgo, 'legend')

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

    # Display the structures next to each other.
    cmd.set('grid_mode')
    # Display the sequences.
    cmd.set('seq_view')
