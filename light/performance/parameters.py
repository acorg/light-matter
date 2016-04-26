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
}
