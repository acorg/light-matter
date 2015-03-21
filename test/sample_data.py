from cStringIO import StringIO

from dark.reads import AARead

from light.database import Database
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.beta_strand import BetaStrand
from light.trig.amino_acids import AminoAcids


DB = Database([AlphaHelix, BetaStrand], [AminoAcids])

_ALPHA = 'ADDDADDDAM'
_BETA = 'VVVVVVM'
_TRYPTOPHAN = 'W'

CATPOX = AARead(
    'Catpox', _ALPHA + _BETA + _ALPHA + _BETA)

COWPOX = AARead(
    'Cowpox', _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN)

MONKEYPOX = AARead(
    'Monkeypox', _BETA + _ALPHA + _ALPHA + _ALPHA + _ALPHA)

MUMMYPOX = AARead(
    'Mummypox', _BETA + _ALPHA + _ALPHA + _ALPHA + _ALPHA + _TRYPTOPHAN)

SQUIRRELPOX = AARead(
    'Squirrelpox', _ALPHA + _BETA + _BETA + _ALPHA + _BETA)

_CATPOX_INDEX = DB.addSubject(CATPOX)
_COWPOX_INDEX = DB.addSubject(COWPOX)
_MONKEYPOX_INDEX = DB.addSubject(MONKEYPOX)
_MUMMYPOX_INDEX = DB.addSubject(MUMMYPOX)
_SQUIRRELPOX_INDEX = DB.addSubject(SQUIRRELPOX)

PARAMS = DB.saveParamsAsJSON(StringIO()).getvalue()

# Run find on a read that matches squirrelpox and catpox.
READ0 = AARead('read0', _ALPHA + _BETA + _BETA + _ALPHA + _BETA)
_result = DB.find(READ0, storeFullAnalysis=True, significanceFraction=0.2)
READ0_SQUIRRELPOX_SCORE = _result.analysis[_SQUIRRELPOX_INDEX]['bestScore']
READ0_CATPOX_SCORE = _result.analysis[_CATPOX_INDEX]['bestScore']
RECORD0 = _result.save(StringIO()).getvalue()

# Run find on a read that matches both monkeypox and mummypox.
READ1 = AARead('read1', _BETA + _ALPHA + _ALPHA + _ALPHA + _BETA + _TRYPTOPHAN)
_result = DB.find(READ1, storeFullAnalysis=True, significanceFraction=0.25)
READ1_MONKEYPOX_SCORE = _result.analysis[_MONKEYPOX_INDEX]['bestScore']
READ1_MONKEYPOX_HSP2_SCORE = _result.analysis[_MONKEYPOX_INDEX][
    'significantBins'][1]['score']
READ1_MUMMYPOX_SCORE = _result.analysis[_MUMMYPOX_INDEX]['bestScore']
RECORD1 = _result.save(StringIO()).getvalue()

# Run find on a read that matches only cowpox.
READ2 = AARead('read2',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_result = DB.find(READ2, storeFullAnalysis=True, significanceFraction=0.3)
READ2_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD2 = _result.save(StringIO()).getvalue()

# Run find on a second read that also matches just cowpox.
READ3 = AARead('read3',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_result = DB.find(READ3, storeFullAnalysis=True, significanceFraction=0.3)
READ3_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD3 = _result.save(StringIO()).getvalue()

# Run find on a third read that also matches just cowpox.
READ4 = AARead('read4',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_result = DB.find(READ4, storeFullAnalysis=True, significanceFraction=0.3)
READ4_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD4 = _result.save(StringIO()).getvalue()
