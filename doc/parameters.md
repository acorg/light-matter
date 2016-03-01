# Parameters

## Database parameters
These parameters are used for the construction of a database.

`landmarks`: The landmarks used to make a database. Either None (to use the default landmark finders) or a mixed list of landmark finder classes or str landmark finder names. To specify no landmark finders, pass an empty list. Default: `AlphaHelix`, `AlphaHelix_3_10`, `AlphaHelix_pi`, `AminoAcids`, `BetaStrand`, `BetaTurn`, `GOR4AlphaHelix`, `GOR4BetaStrand`, `Prosite`

`trigPoints`: The trigPoints used to make a database. Either None (to use the default trig point finders) or a mixed list of trig point finder classes or str trig point finder names. To specify no trig point finders, pass an empty list. Default: `AminoAcids`

`limitPerLandmark`: The maximum number of pairs a landmark can be involved in. Default: 10

`maxDistance`: The maximum distance two features can be apart from each other to form a pair. Default: 200

`minDistance`: The minimum distance two features can be apart from each other to form a pair. Default: 1

`distanceBase`: The pairDistance  between a landmark and a trig point is scaled to be its logarithm using this distanceBase. This reduces sensitivity to relatively small differences in distance. Scaling is done logarithmically to be less sensitive to features which are far apart. Default: 1.1

`featureLengthBase`: The length of a landmark is scaled to be its logarithm using this base, for the purpose of matching landmarks via hashes. This reduces sensitivity to relatively small differences in feature lengths. Default: 1.35

`randomLandmarkDensity`: The density of random length-1 landmarks to generate. This parameter is specific to the RandomLandmark finder. Default: 0.1

`randomTrigPointDensity`: The density of random trig points to generate. This parameter is specific to the RandomTrigPoint finder. Default: 0.1


## Find parameters
These parameters are used for finding and scoring matches.

`significanceMethod`: The name of the method used to calculate which histogram bins are considered significant. This has to be one of `Always`, `HashFraction`, `MaxBinHeight`, `MeanBinHeight`, `AAFraction`. Default: `HashFraction`

`significanceFraction`: The fraction of all pairs for a scanned read that need to fall into the same histogram bucket for that bucket to be considered a significant match with a subject when using the HashFraction significance method. Default: 0.25

`binScoreMethod`: The name of the method used to calculate the score of a histogram bin which is considered significant. This has to be one of `NoneScore`, `MinHashesScore`, `FeatureMatchingScore`, `FeatureAAScore`, `WeightedFeatureAAScore`. Default: `MinHashesScore`

`featureMatchScore`: The contribution (usually positive) to a score when a feature in a query and subject are part of a match when using FeatureMatchingScore as a bin score method. Default: 1.0

`featureMismatchScore`: The contribution (usually negative) to a score when a feature in a query and subject are not part of a match when using FeatureMatchingScore as a bin score method. Default: -1.0

`weights`: A dictionary, where the keys are feature names and the values are weights that should be assigned to that feature in the WeightedFeatureAAScore calculation.

`deltaScale`: The delta between the queryOffset and the subjectOffset of a matching pair is scaled by dividing the delta by the deltaScale to reduce sensitivity. Default: 1.0
