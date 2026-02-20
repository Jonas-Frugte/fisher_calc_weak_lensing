import bispec_cov_solver as ms
import numpy as np

# converting between lists of lists and flat lists

def index_2d_5(X, Y):
    return X * 5 + Y

def index_3d_5(X, Y, Z):
    return X * 25 + Y * 5 + Z

def flatten_2d(matrix):
    matrix_flat = np.zeros(25)
    for i in range(5):
        for j in range(5):
            matrix_flat[index_2d_5(i, j)] = matrix[i, j]
    return matrix_flat

def flatten_3d(matrix):
    matrix_flat = np.zeros(125)
    for i in range(5):
        for j in range(5):
            for k in range(5):
                matrix_flat[index_3d_5(i, j, k)] = matrix[i, j, k]
    return matrix_flat

# testing equation solver

matrix_for_es = np.random.rand(5, 5)

vec = np.random.rand(5)
sol = np.zeros(5)

ms.soe_solver_5_wrap(
    flatten_2d(matrix_for_es),
    vec,
    sol
)
print("\nEquations solver test:")
print("Answer we want: ", vec)
print("Answer we get: ", matrix_for_es @ sol, "\n")

# testing solve_1

CXp = np.random.rand(5, 5)
Yp = np.random.randint(0, 5)
Zp = np.random.randint(0, 5)
bispec_vec = np.random.rand(5, 5, 5)
out5 = np.zeros(5)

ms.solve_1_wrap(
    flatten_2d(CXp),
    Yp,
    Zp,
    flatten_3d(bispec_vec),
    out5
)

print("solve_1 test:")
print("Answer we want: ", np.linalg.solve(CXp, bispec_vec[:, Yp, Zp]))
print("Answer we get: ", out5, "\n")

# testing solve_2

CXp = np.random.rand(5, 5)
CYp = np.random.rand(5, 5)
Zp = np.random.randint(0, 5)
bispec_vec = np.random.rand(5, 5, 5)
out25 = np.zeros(25)

ms.solve_2_wrap(
    flatten_2d(CXp),
    flatten_2d(CYp),
    Zp,
    flatten_3d(bispec_vec),
    out25
)

print("solve_2 test:")
print("Answer we want: ", flatten_2d(np.einsum("ax,by,xyz->abz", np.linalg.inv(CXp), np.linalg.inv(CYp), bispec_vec)[:, :, Zp]))
print("Answer we get: ", out25, "\n")

# testing solve_3

CXp = np.random.rand(5, 5)
CYp = np.random.rand(5, 5)
CZp = np.random.rand(5, 5)
bispec_vec = np.random.rand(5, 5, 5)
out125 = np.zeros(125)

ms.solve_3_wrap(
    flatten_2d(CXp),
    flatten_2d(CYp),
    flatten_2d(CZp),
    flatten_3d(bispec_vec),
    out125
)

print("solve_3 test:")
print("Answer we want: ", flatten_3d(np.einsum("ax,by,cz,xyz->abc", np.linalg.inv(CXp), np.linalg.inv(CYp), np.linalg.inv(CZp), bispec_vec)))
print("Answer we get: ", out125, "\n")