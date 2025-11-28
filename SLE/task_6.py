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

def get_eigenvalues(A, eps, prev_eigenvalues=None):
    n = len(A)
    m = 0
    eigenvalues = []
    end = True
    
    for m in range(n-2):
        sum_squares = 0.0
        for i in range(m + 1, n):
            sum_squares += A[i][m] * A[i][m]

        if sum_squares < eps or m == n-1:
            eigenvalues.append(A[m][m])
            m += 1
        else:
            a, b = A[m][m], A[m][m + 1]
            c, d = A[m + 1][m], A[m + 1][m + 1]
            trace = a + d
            det = a * d - b * c
            disc = trace * trace - 4 * det
            if disc >= 0:
                lambda1 = (trace + math.sqrt(disc)) / 2
                lambda2 = (trace - math.sqrt(disc)) / 2
                eigenvalues.append(lambda1)
                eigenvalues.append(lambda2)
            else:
                real = trace / 2
                imag = math.sqrt(-disc) / 2
                eigenvalues.append(complex(real, imag))
                eigenvalues.append(complex(real, -imag))
            m += 2

            if prev_eigenvalues is None:
                end = False
                continue

            if abs(prev_eigenvalues[m-1] - eigenvalues[m-1]) > eps or abs(prev_eigenvalues[m-1] - eigenvalues[m-1]) > eps:
                end = False
    
    eigenvalues.append(A[n-1][n-1])
    
    return eigenvalues, end

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
        v[0] -= norm_x if x[0] >= 0 else -norm_x

        norm_v = norm(v)
        if norm_v < eps:
            continue

        v = [vi / norm_v for vi in v]

        # H = I - 2vv^T
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
    return eigenvalues

def qr_algorithm(A, max_iterations, eps):
    A_current = copy_matrix(A)
    prev_eigenvalues = None
    eigenvalues = None
    
    for iteration in range(max_iterations):
        Q, R = qr_decomposition(A_current, eps)
        A_new = multiply_matrices(R, Q)

        diag = [A_new[i][i] for i in range(len(A_new))]
        print(f"Итерация {iteration + 1}: диагональ {diag}")

        A_current = A_new
        eigenvalues, stop = get_eigenvalues(A_new, eps, prev_eigenvalues)
        if stop:
            print(f"Сходимость достигнута за {iteration + 1} итераций")
            break

        prev_eigenvalues = eigenvalues
    else:
        print(f"Сходимость достигнута за {max_iterations} итераций")
    print_matrix(A_current)
    return eigenvalues

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

    eps = 0.00001
    eigenvalues = qr_algorithm(A, 1000, eps)

    print("\nСобственные значения:")
    for i, lam in enumerate(eigenvalues, start=1):
        if abs(lam.imag) < eps:
            print(f"λ{i} = {lam.real}")
        else:
            print(f"λ{i} = {lam.real} + {lam.imag}i")

    verify_eigenvalues(A, eigenvalues, eps)
