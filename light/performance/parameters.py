"""
This file holds named parameter sets for use in performance testing.

When running light/performance/bin/perf.py you can refer to these via
--parameterSet XXX (where XXX is a key from the dictionary below).
"""

from light.parameters import DatabaseParameters, FindParameters

PARAMETER_SETS = {
    # The command-line parameter set is special. It is made dynamically
    # from command-line arguments, if any.
    'command-line': {
        'dbParams': None,
        'findParams': None,
    },

    'PDB-landmarks-no-trig': {
        'dbParams': DatabaseParameters(
            landmarks=['PDB AlphaHelix', 'PDB AlphaHelix_3_10',
                       'PDB AlphaHelix_pi', 'PDB ExtendedStrand'],
            trigPoints=[],
            maxDistance=5000,
            limitPerLandmark=50,
        ),
        'findParams': FindParameters(
            binScoreMethod='FeatureAAScore',
            significanceFraction=0.05,
        ),
    },

    'GOR4-landmarks-default-trig': {
        'dbParams': DatabaseParameters(
            landmarks=['GOR4AlphaHelix', 'GOR4BetaStrand', 'GOR4Coil'],
            maxDistance=5000,
            limitPerLandmark=50,
        ),
        'findParams': FindParameters(
            binScoreMethod='FeatureAAScore',
            significanceFraction=0.05,
        ),
    },

    # AC-helices is just the Aho Corasick alpha helix substring finder,
    # with some typical other parameter settings.
    'AC-alpha-helix': {
        'dbParams': DatabaseParameters(
            landmarks=['AC AlphaHelix'],
            trigPoints=[],
            maxDistance=5000,
            limitPerLandmark=50,
            ahocorasickFilename=DatabaseParameters.DEFAULT_AHOCORASICK_FILENAME
        ),
        'findParams': FindParameters(
            binScoreMethod='FeatureAAScore',
            significanceFraction=0.05,
        ),
    },

    'no-finders': {
        # This is included just to allow a quick run of performance tests
        # with very little work being done because no finders are used.
        'dbParams': DatabaseParameters(landmarks=[], trigPoints=[]),
        'findParams': FindParameters(),
    },
}
