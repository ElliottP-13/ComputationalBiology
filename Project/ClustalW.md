# ClustalW Algorithm

1. Calculate all possible pairwise alignments, record the score for each pair. 
2. Calculate a guide tree based on the pairwise distances (algorithm: Neighbor Joining).
3. Find the two most closely related sequences
4. Align the sequences by progressive method
    1. Calculate a consensus of this alignment
    2. Replace the two sequences with the consensus
    3. Find the two next-most closely related sequences (one of these could be a previously determined consensus sequence).
    4. Iterate until all sequences have been aligned

5. Expand the consensus sequences with the (gapped) original sequences
6. Report the multiple sequence alignment