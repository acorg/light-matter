from cStringIO import StringIO

from dark.reads import AARead

from light.database import Database
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi
from light.landmarks.beta_strand import BetaStrand
from light.trig.amino_acids import AminoAcids
from light.trig.troughs import Troughs


DB = Database([AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand],
              [AminoAcids, Troughs])

# An alpha helix with many Tryptophan (W) trig points.
COWPOX = AARead('Cowpox virus 15', 'ADDDADDDAMWDWMWDWMWDW')
DB.addSubject(COWPOX)

# An alpha helix pi with one Tryptophan (W) trig point.
MONKEYPOX = AARead('Monkeypox virus 456', 'ADDDDADDDDAMW')
DB.addSubject(MONKEYPOX)

# A beta strand with one trig point.
MUMMYPOX = AARead('Mummypox virus 3000 B.C.', 'VVVVVVVAW')
DB.addSubject(MUMMYPOX)

# An alpha helix with one trig point.
SQUIRRELPOX1296 = AARead('Squirrelpox virus 1296/99', 'ADDDADDDAMDWMDW')
DB.addSubject(SQUIRRELPOX1296)

# An alpha helix 3, 10 with one trig point.
SQUIRRELPOX55 = AARead('Squirrelpox virus 55', 'ADDADDAMWDW')
DB.addSubject(SQUIRRELPOX55)

PARAMS = DB.saveParamsAsJSON(StringIO()).getvalue()


# The following function is currently unused. It's useful when debugging
# when we change the way we do scoring, significance, etc. For now I'm
# leaving it in as those things are going to change again soon. Eventually
# we can delete it. To use it, add a call to it below after you compute one
# of the _result variables.
def _printResult(result):
    """
    Print information about a database find() result.

    @param result: A C{light.result.Result} instance.
    """
    print 'DB.find() result for read', result.read.id
    for subjectIndex in result.matches:
        print 'subjectIndex', subjectIndex
        print result.analysis[subjectIndex]


# Run find on a read that hits both squirrelpox subjects.
READ0 = AARead('id0', SQUIRRELPOX1296.sequence + SQUIRRELPOX55.sequence)
_result = DB.find(READ0, storeAnalysis=True, significanceFraction=0.18)
RECORD0 = _result.save(StringIO()).getvalue()

# Run find on a read that hits both monkeypox and mummypox. Note that if
# you swap the order of mummypox and monkeypox, you'll need to stick a
# non-beta-sheet AA in the middle (like 'A'), otherwise both subjects will
# not match as the last two residues in the monkeypox sequence will then
# start a beta strand (instead of the mummypox sequence starting it 2
# residues later).
READ1 = AARead('id1', MUMMYPOX.sequence + MONKEYPOX.sequence)
_result = DB.find(READ1, storeAnalysis=True, significanceFraction=0.01)
RECORD1 = _result.save(StringIO()).getvalue()

# Run find on a read that hits only cowpox.
READ2 = AARead('id2', COWPOX.sequence)
_result = DB.find(READ2, storeAnalysis=True, significanceFraction=0.4)
RECORD2 = _result.save(StringIO()).getvalue()

# Run find on a second read that also hits just cowpox.
READ3 = AARead('id3', COWPOX.sequence)
_result = DB.find(READ3, storeAnalysis=True, significanceFraction=0.4)
RECORD3 = _result.save(StringIO()).getvalue()

# Run find on a third read that also hits just cowpox.
READ4 = AARead('id4', COWPOX.sequence)
_result = DB.find(READ4, storeAnalysis=True, significanceFraction=0.4)
RECORD4 = _result.save(StringIO()).getvalue()
