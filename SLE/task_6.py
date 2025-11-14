import math

def print_matrix(A):
    for row in A:
        print(" ".join(f"{x:12.6f}" for x in row))
    print()

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def identity_matrix(n):
    I = [[0.0] * n for _ in range(n)]
    for i in range(n):
        I[i][i] = 1.0
    return I

def multiply_matrices(A, B):
    n, m, p = len(A), len(B), len(B[0])
    C = [[0.0] * p for _ in range(n)]
    for i in range(n):
        for j in range(p):
            for k in range(m):
                C[i][j] += A[i][k] * B[k][j]
    return C

def copy_matrix(A):
    return [row[:] for row in A]

def is_quasi_upper_triangular(A, eps):
    n = len(A)
    for i in range(1, n):
        if i > 0 and abs(A[i][i - 1]) > eps:
            if i == 1 or abs(A[i - 1][i - 2]) < eps:
                continue
        for j in range(0, i - 1):
            if abs(A[i][j]) > eps:
                return False
    return True

def qr_decomposition(A, eps):
    n = len(A)
    Q = identity_matrix(n)
    R = copy_matrix(A)

    for k in range(n - 1):
        # Вектор x из k-го столбца начиная с k-й строки
        x = [R[i][k] for i in range(k, n)]
        norm_x = norm(x)
        if norm_x < eps:
            continue

        v = [xi for xi in x]
        v[0] += norm_x if x[0] >= 0 else -norm_x

        norm_v = norm(v)
        if norm_v < eps:
            continue

        v = [vi / norm_v for vi in v]

        # Матрица Хаусхолдера H = I - 2vv^T
        H = identity_matrix(n)
        for i in range(k, n):
            for j in range(k, n):
                H[i][j] -= 2 * v[i - k] * v[j - k]

        R = multiply_matrices(H, R)
        Q = multiply_matrices(Q, H)
    return Q, R

def determinant(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    det = 0
    for j in range(n):
        minor = [[A[i][k] for k in range(n) if k != j] for i in range(1, n)]
        det += ((-1) ** j) * A[0][j] * determinant(minor)
    return det

def subtract_lambda_from_diagonal(A, lam):
    n = len(A)
    result = copy_matrix(A)
    for i in range(n):
        result[i][i] -= lam
    return result

def verify_eigenvalues(A, eigenvalues, eps):
    print("\nПроверка собственных значений через характеристическое уравнение:")
    print("det(A - λI) должен быть ≈ 0\n")
    for i, lam in enumerate(eigenvalues, start=1):
        if abs(lam.imag) < eps:
            A_minus_lambda_I = subtract_lambda_from_diagonal(A, lam.real)
            det_value = determinant(A_minus_lambda_I)
            print(f"λ{i} = {lam.real:.6f}: det(A - λI) = {det_value:.2e}")
        else:
            print(f"λ{i} = {lam.real:.6f} + {lam.imag:.6f}i: (комплексное, пропущено)")

def extract_eigenvalues(A, eps):
    n = len(A)
    eigenvalues = []
    i = 0
    while i < n:
        if i == n - 1 or abs(A[i + 1][i]) < eps:
            eigenvalues.append(complex(A[i][i], 0))
            i += 1
        else:
            a, b = A[i][i], A[i][i + 1]
            c, d = A[i + 1][i], A[i + 1][i + 1]
            trace = a + d
            det = a * d - b * c
            disc = trace * trace - 4 * det
            if disc >= 0:
                lambda1 = (trace + math.sqrt(disc)) / 2
                lambda2 = (trace - math.sqrt(disc)) / 2
                eigenvalues.extend([complex(lambda1, 0), complex(lambda2, 0)])
            else:
                real = trace / 2
                imag = math.sqrt(-disc) / 2
                eigenvalues.extend([complex(real, imag), complex(real, -imag)])
            i += 2
    return eigenvalues

def qr_algorithm(A, max_iterations, eps):
    A_current = copy_matrix(A)
    for iteration in range(max_iterations):
        Q, R = qr_decomposition(A_current, eps)
        A_new = multiply_matrices(R, Q)

        diag = [A_new[i][i] for i in range(len(A_new))]
        print(f"Итерация {iteration + 1}: диагональ {diag}")

        A_current = A_new
        if is_quasi_upper_triangular(A_new, eps):
            print(f"Сходимость достигнута за {iteration + 1} итераций")
            break
    else:
        print(f"Достигнуто максимальное число итераций: {max_iterations}")
    print_matrix(A_current)
    return extract_eigenvalues(A_current, eps)

if __name__ == "__main__":
    A = [
        [4, -5, -4, 7, 8],
        [4, 10, 0, 4, 2],
        [-2, 3, 4, -10, 5],
        [2, 4, -4, 1, 3],
        [5, 3, -5, 1, -2]
    ]

    print("Исходная матрица:")
    print_matrix(A)

    eps = 0.00000001
    eigenvalues = qr_algorithm(A, 1000, eps)

    print("\nСобственные значения:")
    for i, lam in enumerate(eigenvalues, start=1):
        if abs(lam.imag) < eps:
            print(f"λ{i} = {lam.real}")
        else:
            print(f"λ{i} = {lam.real} + {lam.imag}i")

    verify_eigenvalues(A, eigenvalues, eps)
