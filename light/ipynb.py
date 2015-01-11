# Import and make available things that are convenient to have around in
# iPythonNotebook following 'from light.ipynb import *'.

from dark.fasta import FastaReads
from dark.titles import titleCounts, TitlesAlignments
from dark.html import (
    summarizeTitlesByLength, summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore, summarizeTitlesByCount,
    summarizeTitlesByTitle)
from .alignments import LightReadsAlignments
from .graphics import plotFeatures, plotHistogram

# Keep pyflakes quiet by pretending to make use of all our imports.
_ = (FastaReads, LightReadsAlignments, titleCounts, TitlesAlignments,
     summarizeTitlesByLength, summarizeTitlesByMaxScore,
     summarizeTitlesByMedianScore, summarizeTitlesByCount,
     summarizeTitlesByTitle, plotFeatures)

__all__ = [
    'FastaReads', 'LightReadsAlignments', 'titleCounts', 'TitlesAlignments',
    'summarizeTitlesByLength', 'summarizeTitlesByMaxScore',
    'summarizeTitlesByMedianScore', 'summarizeTitlesByCount',
    'summarizeTitlesByTitle', 'plotFeatures', 'plotHistogram']
