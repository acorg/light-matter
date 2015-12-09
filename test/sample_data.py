from six import StringIO

from dark.reads import AARead

from light.parameters import FindParameters
from light.database import Database, Parameters
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.beta_strand import BetaStrand
from light.trig.amino_acids import AminoAcids


_params = Parameters([AlphaHelix, BetaStrand], [AminoAcids])
DB = Database(_params)

PARAMS = _params.save(StringIO()).getvalue()

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

_, _CATPOX_INDEX, _ = DB.addSubject(CATPOX)
_, _COWPOX_INDEX, _ = DB.addSubject(COWPOX)
_, _MONKEYPOX_INDEX, _ = DB.addSubject(MONKEYPOX)
_, _MUMMYPOX_INDEX, _ = DB.addSubject(MUMMYPOX)
_, _SQUIRRELPOX_INDEX, _ = DB.addSubject(SQUIRRELPOX)

# Run find on a read that matches squirrelpox and catpox.
READ0 = AARead('read0', _ALPHA + _BETA + _BETA + _ALPHA + _BETA)
_findParams = FindParameters(significanceFraction=0.2)
_result = DB.find(READ0, _findParams, storeFullAnalysis=True)
READ0_SQUIRRELPOX_SCORE = _result.analysis[_SQUIRRELPOX_INDEX]['bestScore']
READ0_CATPOX_SCORE = _result.analysis[_CATPOX_INDEX]['bestScore']
RECORD0 = _result.save(StringIO()).getvalue()

# Run find on a read that matches both monkeypox and mummypox.
READ1 = AARead('read1', _BETA + _ALPHA + _ALPHA + _ALPHA + _BETA + _TRYPTOPHAN)
_findParams = FindParameters(significanceFraction=0.25)
_result = DB.find(READ1, _findParams, storeFullAnalysis=True)
READ1_MONKEYPOX_SCORE = _result.analysis[_MONKEYPOX_INDEX]['bestScore']
READ1_MONKEYPOX_HSP2_SCORE = _result.analysis[_MONKEYPOX_INDEX][
    'significantBins'][1]['score']
READ1_MUMMYPOX_SCORE = _result.analysis[_MUMMYPOX_INDEX]['bestScore']
RECORD1 = _result.save(StringIO()).getvalue()

# Run find on a read that matches only cowpox.
READ2 = AARead('read2',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_findParams = FindParameters(significanceFraction=0.3)
_result = DB.find(READ2, _findParams, storeFullAnalysis=True)
READ2_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD2 = _result.save(StringIO()).getvalue()

# Run find on a second read that also matches just cowpox.
READ3 = AARead('read3',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_findParams = FindParameters(significanceFraction=0.3)
_result = DB.find(READ3, _findParams, storeFullAnalysis=True)
READ3_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD3 = _result.save(StringIO()).getvalue()

# Run find on a third read that also matches just cowpox.
READ4 = AARead('read4',
               _ALPHA + _TRYPTOPHAN + _TRYPTOPHAN + _TRYPTOPHAN + _BETA)
_findParams = FindParameters(significanceFraction=0.3)
_result = DB.find(READ4, _findParams, storeFullAnalysis=True)
READ4_COWPOX_SCORE = _result.analysis[_COWPOX_INDEX]['bestScore']
RECORD4 = _result.save(StringIO()).getvalue()
