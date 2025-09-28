import copy

def lu_solve_full(A_in, b_in, eps=1e-12, round_digits=3):
    A = copy.deepcopy(A_in)
    b = copy.deepcopy(b_in)
    n = len(A)

    L = [[0.0]*n for _ in range(n)]
    U = [[0.0]*n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0

    for i in range(n):
        for j in range(i, n):
            s = 0.0
            for k in range(i):
                s += L[i][k] * U[k][j]
            U[i][j] = A[i][j] - s

        for j in range(i+1, n):
            s = 0.0
            for k in range(i):
                s += L[j][k] * U[k][i]
            L[j][i] = (A[j][i] - s) / U[i][i]

    det = 1.0
    for i in range(n):
        det *= U[i][i]

    y = [0.0]*n
    for i in range(n):
        s = 0.0
        for k in range(i):
            s += L[i][k] * y[k]
        y[i] = b[i] - s

    x = [0.0]*n
    for i in range(n-1, -1, -1):
        s = 0.0
        for k in range(i+1, n):
            s += U[i][k] * x[k]
        x[i] = (y[i] - s) / U[i][i]

    inv = [[0.0]*n for _ in range(n)]

    for col in range(n):
        e = [0.0]*n
        e[col] = 1.0
        ycol = [0.0]*n
        for i in range(n):
            s = 0.0
            for k in range(i):
                s += L[i][k] * ycol[k]
            ycol[i] = e[i] - s
        xcol = [0.0]*n
        for i in range(n-1, -1, -1):
            s = 0.0
            for k in range(i+1, n):
                s += U[i][k] * xcol[k]
            xcol[i] = (ycol[i] - s) / U[i][i]
        for i in range(n):
            inv[i][col] = xcol[i]

    x_rounded = [round(v, round_digits) for v in x]
    det_rounded = round(det, round_digits)
    inv_rounded = [[round(v, round_digits) for v in row] for row in inv]

    return {
        "x": x, "det": det, "inv": inv,
        "x_rounded": x_rounded, "det_rounded": det_rounded, "inv_rounded": inv_rounded,
        "L": L, "U": U
    }


if __name__ == "__main__":
    A = [
        [4, -6, -3, 9, 4],
        [-5, 3, 1, 4, 2],
        [-2, 1, 4, 2, -5],
        [7, 5, 0, 1, 3],
        [8, -3, 4, 1, -2]
    ]
    b = [31, 62, 47, 77, 52]

    res = lu_solve_full(A, b, round_digits=3)

    print("Решение x:", res["x_rounded"])
    print("Детерминант:", res["det_rounded"])
    print("Обратная матрица:")
    for row in res["inv_rounded"]:
        print(row)

