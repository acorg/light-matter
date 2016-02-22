# Definitions of terms

`Pair`: Two features which are considered together. Pairs can consist of either two landmarks or one landmark and one trigPoint.

`pairDistance`: The distance between the two features in a pair.

`delta`: When two pairs match, the delta is the difference between the distance from the start of the subject to the offset of the pair in the subject and the distance from the start of the query to the offset of the pair in the query.

`feature offset`: The offset of a feature in the sequence.

`bin`: All deltas from the comparison of two sequences are entered in a histogram. A bin thus contans deltas of similar length.

`matching region`: The matching region is the region from the smallest feature offset in the bin to the largest feature offset in the bin.
