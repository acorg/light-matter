# Sample light matter parameters.
PARAMS = {
    'checksum': 'efb324c955afe55b485efcb144e9592911c014eb7fd77acd4d5bc672d68b',
    'landmarkFinderClasses': [
        'AlphaHelix_pi',
        'AlphaHelix',
        'AlphaHelix_3_10',
        'BetaStrand'
    ],
    'limitPerLandmark': 10,
    'maxDistance': 50,
    'subjectCount': 20151,
    'totalCoveredResidues': 949153,
    'totalResidues': 1252998,
    'trigPointFinderClasses': [
        'AminoAcids',
        'IndividualPeaks',
        'IndividualTroughs',
        'Troughs',
        'Peaks'
    ]
}

RECORD0 = {
    'alignments': [
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 60,
            'subjectIndex': 0,
        },
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 50,
            'subjectIndex': 1,
        }
    ],
    'query': 'SK7F6:446:2043-frame0rc-[0:60)',
    'querySequence': 'ARLRYVQTCAMHPSMNFPRLSHVRSLTYSDKPRAQHFYFDNP'
}

RECORD1 = {
    'alignments': [
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 60,
            'subjectIndex': 0,
        },
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 50,
            'subjectIndex': 1,
        }
    ],
    'query': 'SK7F6:446:2043-frame0rc-[0:60)',
    'querySequence': 'ARLRYVQTCAMHPSMNFPRLSHVRSLTYSDKPRAQHFYFDNP'
}

RECORD2 = {
    'alignments': [
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 60,
            'subjectIndex': 0,
        },
    ],
    'query': 'SK7F6:446:2043-frame0rc-[0:60)',
    'querySequence': 'ARLRYVQTCAMHPSMNFPRLSHVRSLTYSDKPRAQHFYFDNP'
}

# Identical to RECORD2, apart from match score.
RECORD3 = {
    'alignments': [
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
            ],
            'matchScore': 70,
            'subjectIndex': 0,
        },
    ],
    'query': 'SK7F6:446:2043-frame0rc-[0:60)',
    'querySequence': 'ARLRYVQTCAMHPSMNFPRLSHVRSLTYSDKPRAQHFYFDNP'
}

RECORD4 = {
    'alignments': [
        {
            'hsps': [
                {
                    'readOffset': 27,
                    'subjectOffset': 27
                },
                {
                    'readOffset': 27,
                    'subjectOffset': 43
                },
                {
                    'readOffset': 37,
                    'subjectOffset': 48
                },
            ],
            'matchScore': 70,
            'subjectIndex': 0,
        },
    ],
    'query': 'SK7F6:446:2043-frame0rc-[0:60)',
    'querySequence': 'ARLRYVQTCAMHPSMNFPRLSHVRSLTYSDKPRAQHFYFDNP'
}
