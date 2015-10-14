#!/usr/bin/env python

import sys
import argparse
from time import time
from json import loads
from os.path import dirname, join
import matplotlib.pyplot as plt
from scipy import stats

from dark.fasta import FastaReads
from dark.reads import AARead

import light
from light.database import DatabaseSpecifier


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

DATASETS = {
    '2VQ0': VQ,
    '3NXG': NXG,
    '4MTP': MTP,
}


def makeScatterplot(lightScores, otherScores, scoreType, outFile, dataset,
                    regression=False):
    """
    Make a scatterplot of the light matter score and another score.

    @param lightScores: A C{list} of light matter scores.
    @param otherScores: A C{list} of other scores.
    @param scoreType: A C{str} describing the type of scores in scores. Must be
        one of 'pid', 'z', 'rmsd'.
    @param outFile: The C{str} name of the file where the figure will be
        written to.
    @param dataset: The C{str} name of the dataset that was used to make the
        plot.
    @param regression: If True, a regression line will be added to the figure.
    """
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    plt.plot(lightScores, otherScores, 'o', markerfacecolor='blue',
             markeredgecolor='white')

    if regression:
        slope, intercept, rValue, pValue, se = stats.linregress(lightScores,
                                                                otherScores)
        plt.plot([0, max(lightScores)], [intercept,
                 slope * max(lightScores) + intercept], '-', color='black')

    # labels
    ax.set_xlabel('Light matter score')
    ax.set_ylabel(scoreType)
    if regression:
        ax.set_title('%s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' % (
                     dataset, rValue, se, slope, pValue))
    else:
        ax.set_title(dataset)

    # axes
    if scoreType == 'z':
        ax.set_ylim(0, 70)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    fig.savefig(outFile)
    print('Wrote scatterplot to %s.' % outFile, file=sys.stderr)


if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description='Create a scatter plot of the light matter score and '
        'either percent sequence identity, Z-score or RMSD.')

    parser.add_argument(
        '--dataset', choices=(DATASETS.keys()), required=True,
        help='The name of the dataset which should be used.')

    parser.add_argument(
        '--outFile', required=True,
        help='The name of the output file.')

    parser.add_argument(
        '--score', choices=('pid', 'rmsd', 'z'), required=True,
        help='The name of the dataset which should be used.')

    parser.add_argument(
        '--regression', default=False, action='store_true',
        help='If True, add a regression line to the plot.')

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)

    args = parser.parse_args()

    # Get the names of the data files.
    databaseFile = join(dirname(light.__file__), '..', 'data',
                        '%s-db.fasta' % args.dataset)
    daliResultFile = join(dirname(light.__file__), '..', 'data',
                          '%s.json' % args.dataset)

    read = DATASETS.get(args.dataset.upper())
    if read is None:
        raise KeyError('%s is not a valid dataset.' % args.dataset)

    database = databaseSpecifier.getDatabaseFromArgs(args)
    list(map(database.addSubject, FastaReads(databaseFile, readClass=AARead,
                                             upperCase=True)))

    lookupTime = time()
    result = database.find(read)
    print('Look up done in %.2f seconds.' % (time() - lookupTime),
          file=sys.stderr)

    # Read the daliResultFile.
    allScores = {}
    with open(daliResultFile) as fp:
        for line in fp:
            scores = loads(line)
            allScores[scores['title']] = {'z': scores['z'],
                                          'rmsd': scores['rmsd'],
                                          'pid': scores['pid'],
                                          'lm': 0}
    for subjectIndex in result.significantSubjects():
        lightScore = result.analysis[subjectIndex]['score']
        lightTitle = database.getSubjectByIndex(subjectIndex).id
        allScores[lightTitle]['lm'] = lightScore

    # Collect data for plotting.
    lightScores = []
    otherScores = []
    for title in allScores:
        lightScores.append(allScores[title]['lm'])
        otherScores.append(allScores[title][args.score])

    makeScatterplot(lightScores, otherScores, args.score, args.outFile,
                    args.dataset, regression=args.regression)
