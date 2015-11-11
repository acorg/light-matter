import builtins
from io import StringIO
from unittest import TestCase
from unittest.mock import patch
from .mocking import mockOpen
from json import load

from light.prosite import prositeToJSON, patternToRegex


class TestProsite(TestCase):
    """
    Tests for the light.prosite.prositeToJSON and light.prosite.patternToRegex
    functions.
    """
    def testRecordOutput(self):
        """
        A valid prosite record must result in the desired JSON.
        """
        data = '\n'.join(['//',
                          'ID   ASN_GLYCOSYLATION; PATTERN.',
                          'AC   PS00001;',
                          'DT   APR-1990 (CREATED); APR-1990 (DATA UPDATE); '
                          'APR-1990 (INFO UPDATE).',
                          'DE   N-glycosylation site.',
                          'PA   N(3,4)-{P}-[ST]-{P}.',
                          'CC   /SITE=1,carbohydrate;',
                          'CC   /SKIP-FLAG=TRUE;',
                          'CC   /VERSION=1;'
                          'PR   PRU00498;',
                          'DO   PDOC00001;',
                          '//'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            fp = StringIO()
            prositeToJSON('prosite.dat', fp=fp)
            fp.seek(0)
            result = load(fp)
            self.assertEqual({'pattern': 'N{3,4}[^P][ST][^P]',
                              'accession': '00001'},
                             result)

    def testPatternToRegex(self):
        """
        The prosite pattern must be translated to the right regex.
        """
        pattern = ('D-{W}-[DNS](3,7)-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-'
                   '[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]')
        regex = patternToRegex(pattern)
        expected = ('D[^W][DNS]{3,7}[^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC]'
                    '[DENQSTAGC].{2}[DE][LIVMFYW]')
        self.assertEqual(expected, regex)
