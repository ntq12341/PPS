import numpy as np

with open("inp", "r") as f:
    n = int(f.readline())
    A = np.array([list(map(float, f.readline().split())) for _ in range(n)])

eigenvalues, eigenvectors = np.linalg.eig(A)

with open("out", "w") as f:
    f.write("Gia tri rieng:\n")
    for val in eigenvalues:
        f.write(f"{val:.6f} ")
    f.write("\n\nVector rieng:\n")

    for i in range(n):
        for j in range(n):
            f.write(f"{eigenvectors[i, j]:.6f} ")
        f.write("\n")
