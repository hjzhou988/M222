import numpy as np

# alignment matrix:
def alignment_matrix(ref,donor):
    import numpy as np
    output_matrix = np.zeros((len(ref), len(donor)),dtype=int)
    mismatch=-1
    match=1
    insertion=-3
    deletion=-3
    for i in range(len(ref)):
        output_matrix[i, 0] = deletion*i
    for j in range(len(donor)):
        output_matrix[0, j] = deletion*j
    for j in range(1, len(donor)):
        for i in range(1, len(ref)):  # Big opportunities for improvement right here.
            deletion = output_matrix[i - 1, j] + deletion #1. Make it bigger or smaller?
            insertion = output_matrix[i, j - 1] + insertion # 1
            identity = output_matrix[i - 1, j - 1] + match  if ref[i] == donor[j] else -np.inf
            substitution = output_matrix[i - 1, j - 1] + mismatch if ref[i] != donor[j] else -np.inf
            output_matrix[i, j] = max(insertion, deletion, identity, substitution)
    return output_matrix   

alignment_matrix('$TAACGCCCATC','$GCCC')
print('hello world')