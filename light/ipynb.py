# Import and make available things that are convenient to have around in
# iPythonNotebook following 'from light.ipynb import *'.

from dark.reads import AARead
from dark.fasta import FastaReads
from dark.titles import titleCounts, TitlesAlignments
from dark.html import (
    summarizeTitlesByLength, summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore, summarizeTitlesByCount,
    summarizeTitlesByTitle)

from .alignments import LightReadsAlignments
from .database import Database, DatabaseSpecifier
from .graphics import (plotHistogram, plotFeatureSquare)

# Keep pyflakes quiet by pretending to make use of all our imports.
_ = (AARead, FastaReads, LightReadsAlignments, titleCounts, TitlesAlignments,
     summarizeTitlesByLength, summarizeTitlesByMaxScore,
     summarizeTitlesByMedianScore, summarizeTitlesByCount,
     summarizeTitlesByTitle, plotHistogram, plotFeatureSquare, Database,
     DatabaseSpecifier)

__all__ = [
    'AARead', 'FastaReads', 'LightReadsAlignments', 'titleCounts',
    'TitlesAlignments', 'summarizeTitlesByLength', 'summarizeTitlesByMaxScore',
    'summarizeTitlesByMedianScore', 'summarizeTitlesByCount',
    'summarizeTitlesByTitle', 'plotHistogram', 'plotFeatureSquare', 'Database',
    'DatabaseSpecifier']
