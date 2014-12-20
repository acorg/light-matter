# from pprint import pprint
from cStringIO import StringIO

from dark.reads import AARead

from light.database import Database
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi
from light.landmarks.beta_strand import BetaStrand
from light.trig.amino_acids import AminoAcids


DB = Database([AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand],
              [AminoAcids])

# An alpha helix with one (distant) trig point.
COWPOX = AARead('Cowpox virus 15', 'FRRRFRRRFAAAAAAAAAAAAAAAAAAAAAAC')

# An alpha helix pi with one trig point.
MONKEYPOX = AARead('Monkeypox virus 456', 'FRRRRFRRRRFAC')

# A beta strand with one trig point.
MUMMYPOX = AARead('Mummypox virus 3000 B.C.', 'VVVVVVAC')

# An alpha helix with one trig point.
SQUIRRELPOX1296 = AARead('Squirrelpox virus 1296/99', 'FRRRFRRRFAC')

# An alpha helix 3, 10 with one trig point.
SQUIRRELPOX55 = AARead('Squirrelpox virus 55', 'FRRFRRFAC')

DB.addSubject(COWPOX)
DB.addSubject(MONKEYPOX)
DB.addSubject(MUMMYPOX)
DB.addSubject(SQUIRRELPOX1296)
DB.addSubject(SQUIRRELPOX55)

PARAMS = DB.saveParamsAsJSON(StringIO()).getvalue()

# Run find on a read that hits both squirrelpox subjects.
READ0 = AARead('id0', SQUIRRELPOX1296.sequence + SQUIRRELPOX55.sequence)
_result = DB.find(READ0, aboveMeanThreshold=1)
# print 'find result for', SQUIRRELPOX1296.id
# for subjectIndex in _result.matches:
#     print '>>> SUBJECT: %r (index %d)' % (DB.subjectInfo[subjectIndex][0],
#                                           subjectIndex)
#     pprint(_result.matches[subjectIndex])
# print 'significant:'
# pprint(_result.significant)
RECORD0 = _result.save(StringIO()).getvalue()

# Run find on a read that hits both monkeypox and mummypox.
READ1 = AARead('id1', MONKEYPOX.sequence + MUMMYPOX.sequence)
_result = DB.find(READ1, aboveMeanThreshold=1)
RECORD1 = _result.save(StringIO()).getvalue()

# Run find on a read that hits just cowpox.
READ2 = AARead('id2', COWPOX.sequence)
_result = DB.find(READ2, aboveMeanThreshold=1)
RECORD2 = _result.save(StringIO()).getvalue()

# Run find on a second read that also hits just cowpox.
READ3 = AARead('id3', COWPOX.sequence)
_result = DB.find(READ3, aboveMeanThreshold=1)
RECORD3 = _result.save(StringIO()).getvalue()

# Run find on a third read that also hits just cowpox.
READ4 = AARead('id4', COWPOX.sequence)
_result = DB.find(READ4, aboveMeanThreshold=1)
RECORD4 = _result.save(StringIO()).getvalue()
