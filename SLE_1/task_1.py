import copy

def gauss_full(A, b):
    n = len(A)

    M = [A[i] + [b[i]] for i in range(n)]
    M_inv = [A[i] + [int(i == j) for j in range(n)] for i in range(n)]

    det = 1 

    for k in range(n):
        pivot = M[k][k]
        det *= pivot

        for j in range(k, n+1):
            M[k][j] /= pivot
        for j in range(2*n):
            M_inv[k][j] /= pivot

        for i in range(k+1, n):
            factor = M[i][k]
            for j in range(k, n+1):
                M[i][j] -= factor * M[k][j]
            for j in range(2*n):
                M_inv[i][j] -= factor * M_inv[k][j]

    x = [0]*n
    for i in range(n-1, -1, -1):
        s = sum(M[i][j] * x[j] for j in range(i+1, n))
        x[i] = M[i][n] - s

    for k in range(n-1, -1, -1):
        for i in range(k-1, -1, -1):
            factor = M_inv[i][k]
            for j in range(2*n):
                M_inv[i][j] -= factor * M_inv[k][j]

    inv = [row[n:] for row in M_inv]

    return x, det, inv


A = [
    [4, -6, -3, 9, 4],
    [-5, 3, 1, 4, 2],
    [-2, 1, 4, 2, -5],
    [7, 5, 0, 1, 3],
    [8, -3, 4, 1, -2]
]
b = [31, 62, 47, 77, 52]

x, detA, invA = gauss_full(A, b)

print("Решение системы:")
print([round(val, 3) for val in x])

print("\nОпределитель:")
print(round(detA, 3))

print("\nОбратная матрица:")
for row in invA:
    print([round(val, 3) for val in row])