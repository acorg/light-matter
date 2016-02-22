# Import and make available things that are convenient to have around in
# iPythonNotebook following 'from light.ipynb import *'.

from dark.reads import AARead, AAReadWithX, SSAARead
from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.titles import titleCounts, TitlesAlignments
from dark.html import (
    summarizeTitlesByLength, summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore, summarizeTitlesByCount,
    summarizeTitlesByTitle)

from .alignments import LightReadsAlignments
from .database import Database, DatabaseSpecifier
from .graphics import (
    plotHistogramPanel, plotHistogram, plotHistogramLinePanel,
    plotHistogramLine, plotHistogramLines, plotFeatureSquare,
    plotHorizontalPairPanel, PlotHashesInSubjectAndRead,
    plotLandmarksInSequences, confusionMatrix, featureComparison,
    SequenceFeatureAnalysis, compareScores, scoreHeatmap,
    alignmentGraph, alignmentGraphMultipleQueries, alignmentPanel)
from .parameters import FindParameters, DatabaseParameters

# Keep pyflakes quiet by pretending to make use of all our imports.
_ = (
    # From dark.reads.
    AARead, AAReadWithX, SSAARead,

    # From dark.fasta.
    FastaReads,

    # From dark.fasta_ss.
    SSFastaReads,

    # From dark.titles.
    titleCounts, TitlesAlignments,

    # From html.
    summarizeTitlesByLength, summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore, summarizeTitlesByCount,

    # From alignments.
    LightReadsAlignments,

    # From database.
    summarizeTitlesByTitle, Database, DatabaseSpecifier,

    # From graphics.
    plotHistogramPanel, plotHistogram, plotHistogramLinePanel,
    plotHistogramLine, plotHistogramLines, plotFeatureSquare,
    plotHorizontalPairPanel, PlotHashesInSubjectAndRead,
    plotLandmarksInSequences, confusionMatrix, featureComparison,
    SequenceFeatureAnalysis, compareScores, scoreHeatmap,
    alignmentGraph, alignmentGraphMultipleQueries, alignmentPanel,

    # From parameters.
    FindParameters, DatabaseParameters,
)

__all__ = [
    # From dark.reads.
    'AARead', 'AAReadWithX',

    # From dark.fasta.
    'FastaReads',

    # From dark.fasta_ss.
    'SSFastaReads',

    # From dark.titles.
    'titleCounts', 'TitlesAlignments',

    # From html.
    'summarizeTitlesByLength', 'summarizeTitlesByMaxScore',
    'summarizeTitlesByMedianScore', 'summarizeTitlesByCount',
    'summarizeTitlesByTitle',

    # From alignments.
    'LightReadsAlignments',

    # From database.
    'Database', 'DatabaseSpecifier',

    # From graphics.
    'plotHistogramPanel', 'plotHistogram', 'plotHistogramLinePanel',
    'plotHistogramLine', 'plotHistogramLines', 'plotFeatureSquare',
    'plotHorizontalPairPanel', 'PlotHashesInSubjectAndRead',
    'plotLandmarksInSequences', 'confusionMatrix', 'featureComparison',
    'SequenceFeatureAnalysis', 'compareScores', 'scoreHeatmap',
    'alignmentGraph', 'alignmentGraphMultipleQueries', 'alignmentPanel',
]
