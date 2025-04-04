import itertools
import numpy as np

mat = np.zeros((3, 3))

for i, j in itertools.product(range(3), repeat = 2):
    if i >= j:
        x = np.random.rand()
        mat[i, j] = x
        # mat[j, i] = x
print(mat)