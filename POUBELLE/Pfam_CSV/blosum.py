"""
Get the official BLOSUM 62 matrix size for the 20 standard amino-acids
"""

# Code source: Patrick Kunzmann
# License: BSD 3 clause

# import numpy as np
# import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
# import biotite.sequence.phylo as phylo
# import biotite.sequence.graphics as graphics

# Obtain BLOSUM62
matrix = align.SubstitutionMatrix.std_protein_matrix()
# print(matrix)


# Matrix should not contain ambiguous symbols or stop signal
matrix = align.SubstitutionMatrix(
    seq.Alphabet(matrix.get_alphabet1().get_symbols()[:-4]),
    seq.Alphabet(matrix.get_alphabet2().get_symbols()[:-4]),
    matrix.score_matrix()[:-4, :-4]
)
similarities = matrix.score_matrix()
print(matrix)
