#!/usr/bin/env python

import sys
import argparse
from time import time
from json import loads
from os.path import dirname, join, basename
import matplotlib.pyplot as plt
from scipy import stats

from dark.fasta import FastaReads
from dark.reads import AARead

import light
from light.landmarks import findLandmark, DEFAULT_LANDMARK_FINDER_CLASSES
from light.trig import findTrigPoint, DEFAULT_TRIG_FINDER_CLASSES
from light.database import Database


# Data
VQ = AARead(("gi|188036137|pdb|2VQ0_A|Chain A, Capsid Structure Of Sesbania "
             "Mosaic Virus Coat Protein Deletion Mutant Rcp(Delta 48 To 59)., "
             "['Viruses', 'ssRNA positive-strand viruses, no DNA stage', "
             "'Sobemovirus']"),
            ("MAKRLSKQQLAKAIANTLETPPQPKAGRRRNRRRQRSAVQQLQPTQAVRIRNPAVSSSRGGITV"
             "LTHSELSAEIGVTDSIVVSSELVMPYTVGTWLRGVAANWSKYSWLSVRYTYIPSCPSSTAGSIH"
             "MGFQYDMADTVPVSVNQLSNLRGYVSGQVWSGSAGLCFINGTRCSDTSTAISTTLDVSKLGKKW"
             "YPYKTSADYATAVGVDVNIATPLVPARLVIALLDGSSSTAVAAGRIYCTYTIQMIEPTASALNN"
             ))

NXG = AARead(("gi|312207952|pdb|3NXG_A|Chain A, Jc Polyomavirus Vp1., "
              "['Viruses', 'dsDNA viruses, no RNA stage', 'Polyomaviridae', "
              "'Polyomavirus']"),
             ("GSHMGGVEVLEVKTGVDSITEVECFLTPEMGDPDEHLRGFSKSISISDTFESDSPNRDMLPCY"
              "SVARIPLPNLNEDLTCGNILMWEAVTLKTEVIGVTSLMNVHSNGQATHDNGAGKPVQGTSFHF"
              "FSVGGEALELQGVLFNYRTKYPDGTIFPKNATVQSQVMNTEHKAYLDKNKAYPVECWVPDPTR"
              "NENTRYFGTLTGGENVPPVLHITNTATTVLLDEFGVGPLCKGDNLYLSAVDVCGMFTNRSGSQ"
              "QWRGLSRYFKVQLRKRRVKN"))

MTP = AARead(("gi|568786774|pdb|4MTP_A|Chain A, Rdrp From Japanesese "
              "Encephalitis Virus., ['Viruses', 'ssRNA positive-strand "
              "viruses, no DNA stage', 'Flaviviridae', 'Flavivirus', "
              "'Japanese encephalitis virus group']"),
             ("VHSNQEKIKKRIQKLKEEFATTWHKDPEHPYRTWTYHGSYEVKATGSASSLVNGVVKLMSKPW"
              "DAIANVTTMAMTDTTPFGQQRVFKEKVDTKAPEPPAGVREVLNETTNWLWAHLSREKRPRLCT"
              "KEEFIKKVNSNAALGAVFAEQNQWSTAREAVNDPRFWEMVDEERENHLRGECHTCIYNMMGKR"
              "EKKPGEFGKAKGSRAIWFMWLGARYLEFEALGFLNEDHWLSRENSGGGVEGSGVQKLGYILRD"
              "IAGKQGGKMYADDTAGWDTRITRTDLENEAKVLELLDGEHRMLARAIIELTYRHKVVKVMRPA"
              "AEGKTVMDVISREDQRGSGQVVTYALNTFTNIAVQLVRLMEAEGVIGPQHLEQLPRKNKIAVR"
              "TWLFENGEERVTRMAISGDDCVVKPLDDRFATALHFLNAMSKVRKDIQEWKPSHGWHDWQQVP"
              "FCSNHFQEIVMKDGRSIVVPCRGQDELIGRARISPGAGWNVKDTACLAKAYAQMWLLLYFHRR"
              "DLRLMANAICSAVPVDWVPTGRTSWSIHSKGEWMTTEDMLQVWNRVWIEENEWMMDKTPIASW"
              "TDVPYVGKREDIWCGSLIGTRSRATWAENIYAAINQVRAVIGKENYVDYMTSLRRYEDVLIQE"
              "DRVI"))


def makeScatterplot(lmScores, scores, scoreType, outFile, regression=False):
    """
    Make a scatterplot of the light matter score and another score.

    @param lmScores: A C{list} of light matter scores.
    @param scores: A C{list} of other scores.
    @param scoreType: A C{str} describing the type of scores in scores. Must be
        one of 'pid', 'z', 'rmsd'.
    @param outFile: The C{str} name of the file where the figure will be
        written to.
    @param regression: If True, a regression line will be added to the figure.
    """
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    plt.plot(lm, s, 'o', markerfacecolor='blue', markeredgecolor='white')

    if args.regression:
        slope, intercept, rValue, pValue, se = stats.linregress(lm, s)
        plt.plot([0, max(lm)], [intercept, slope * max(lm) + intercept], '-',
                 color='black')

    # labels
    ax.set_xlabel('Light matter score')
    ax.set_ylabel(args.score)
    if args.regression:
        ax.set_title('%s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' % (
                     args.dataset, rValue, se, slope, pValue))
    else:
        ax.set_title(args.dataset)

    # axes
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    fig.savefig(args.outFile)
    print >>sys.stderr, 'Wrote scatterplot to %s.' % args.outFile

if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description='Create a scatter plot of the light matter score and '
        'either percent sequence identity, Z-score or RMSD.')

    parser.add_argument(
        '--dataset', choices=('2VQ0', '3NXG', '4MTP'),
        help='The name of the dataset which should be used.')

    parser.add_argument(
        '--outFile',
        help='The name of the output file.')

    parser.add_argument(
        '--landmark', action='append', dest='landmarkFinderNames',
        choices=sorted(cl.NAME for cl in DEFAULT_LANDMARK_FINDER_CLASSES),
        help='The name of a landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', dest='trigFinderNames',
        choices=sorted(cl.NAME for cl in DEFAULT_TRIG_FINDER_CLASSES),
        help='The name of a trig point finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--limitPerLandmark', type=int, default=None,
        help='A limit on the number of pairs to yield per landmark per read.')

    parser.add_argument(
        '--maxDistance', type=int, default=None,
        help='The maximum distance permitted between yielded pairs.')

    parser.add_argument(
        '--minDistance', type=int, default=None,
        help='The minimum distance permitted between yielded pairs.')

    parser.add_argument(
        '--bucketFactor', type=int, default=1,
        help=('A factor by which the distance between landmark and trig point '
              'is divided.'))

    parser.add_argument(
        '--score', choices=('pid', 'rmsd', 'z'),
        help='The name of the dataset which should be used.')

    parser.add_argument(
        '--regression', default=False, action='store_true',
        help='If True, add a regression line to the plot.')

    args = parser.parse_args()

    landmarkFinderNames = (args.landmarkFinderNames or
                           [klass.NAME for klass in
                            DEFAULT_LANDMARK_FINDER_CLASSES])
    trigFinderNames = (args.trigFinderNames or
                       [klass.NAME for klass in DEFAULT_TRIG_FINDER_CLASSES])

    if len(landmarkFinderNames) + len(trigFinderNames) == 0:
        print >>sys.stderr, ('You must specify either landmark or trig point '
                             'finders to find.\n%s') % parser.format_usage()
        sys.exit(1)

    # Make sure all landmark finders requested exist.
    landmarkFinderClasses = []
    for landmarkFinderName in landmarkFinderNames:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinderClasses.append(landmarkFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Make sure all trig point finders requested exist.
    trigFinderClasses = []
    for trigFinderName in trigFinderNames:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinderClasses.append(trigFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), trigFinderName)
            sys.exit(1)

    # Get the names of the data files
    databaseFile = join(dirname(light.__file__), '..', 'data',
                        '%s-db.fasta' % args.dataset)
    daliResultFile = join(dirname(light.__file__), '..', 'data',
                          '%s.json' % args.dataset)

    if args.dataset == '2VQ0':
        read = VQ
    elif args.dataset == '3NXG':
        read = NXG
    elif args.dataset == '4MTP':
        read = MTP

    startTime = time()
    # Create the database, add reads to it.
    database = Database(landmarkFinderClasses, trigFinderClasses,
                        args.limitPerLandmark, args.maxDistance,
                        args.minDistance, args.bucketFactor)

    reads = FastaReads(databaseFile)

    for r in reads:
        database.addSubject(r)

    print >>sys.stderr, 'Database built in %.2f seconds.' % (time() -
                                                             startTime)
    lookupTime = time()
    result = database.find(read)
    print >>sys.stderr, 'Look up done in %.2f seconds.' % (time() - lookupTime)

    # read the daliResultFile
    allScores = {}
    with open(daliResultFile) as fp:
        for line in fp:
            scores = loads(line)
            allScores[scores['title']] = {'z': scores['z'],
                                          'rmsd': scores['rmsd'],
                                          'pid': scores['pid'],
                                          'lm': 0}
    for subjectIndex in result.significant():
        lmScore = result.analysis[subjectIndex]['score']
        lmTitle = database.subjectInfo[subjectIndex][0]
        allScores[lmTitle]['lm'] = lmScore

    # data for plotting
    lm = []
    s = []
    for t in allScores:
        lm.append(allScores[t]['lm'])
        s.append(allScores[t][args.score])

    makeScatterplot(lm, s, args.score, args.outFile,
                    regression=args.regression)
