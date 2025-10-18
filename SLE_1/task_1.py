import copy

def gauss_full(A, b):
    n = len(A)

    print("=" * 60)
    print("МЕТОД ГАУССА")
    print("=" * 60)
    
    print(f"\nРазмер системы: {n}×{n}")
    print("\nИсходная матрица A:")
    for row in A:
        print([f"{x:6.2f}" for x in row])
    print(f"\nВектор b: {b}")

    M = [A[i] + [b[i]] for i in range(n)]
    M_inv = [A[i] + [int(i == j) for j in range(n)] for i in range(n)]

    det = 1 

    print("\n" + "=" * 60)
    print("ПРЯМОЙ ХОД")
    print("=" * 60)

    for k in range(n):
        print(f"\n--- Шаг {k+1} ---")
        pivot = M[k][k]
        det *= pivot

        print(f"Текущий определитель: {det:.6f}")

        for j in range(k, n+1):
            M[k][j] /= pivot
        for j in range(2*n):
            M_inv[k][j] /= pivot

        print(f"После нормировки строки {k}:")
        print("M = [A|b]:")
        for i in range(n):
            print(f"  {i}: {[f'{M[i][j]:8.4f}' for j in range(n+1)]}")

        for i in range(k+1, n):
            factor = M[i][k]
            print(f"Исключение из строки {i} (множитель: {factor:.6f})")
            
            for j in range(k, n+1):
                M[i][j] -= factor * M[k][j]
            for j in range(2*n):
                M_inv[i][j] -= factor * M_inv[k][j]

        print(f"После исключения:")
        print("M = [A|b]:")
        for i in range(n):
            print(f"  {i}: {[f'{M[i][j]:8.4f}' for j in range(n+1)]}")

    print("\n" + "=" * 60)
    print("ОБРАТНЫЙ ХОД - РЕШЕНИЕ")
    print("=" * 60)

    x = [0]*n
    for i in range(n-1, -1, -1):
        s = sum(M[i][j] * x[j] for j in range(i+1, n))
        x[i] = M[i][n] - s
        print(f"x[{i}] = {M[i][n]:.6f} - {s:.6f} = {x[i]:.6f}")

    print("\n" + "=" * 60)
    print("ОБРАТНЫЙ ХОД - ОБРАТНАЯ МАТРИЦА")
    print("=" * 60)

    for k in range(n-1, -1, -1):
        print(f"\n--- Обратный шаг для строки {k} ---")
        for i in range(k-1, -1, -1):
            factor = M_inv[i][k]
            print(f"Строка {i}: вычитаем {factor:.6f} × строку {k}")
            for j in range(2*n):
                M_inv[i][j] -= factor * M_inv[k][j]

        print("M_inv после шага:")
        for i in range(n):
            print(f"  {i}: {[f'{M_inv[i][j]:8.4f}' for j in range(2*n)]}")

    inv = [row[n:] for row in M_inv]

    print("\n" + "=" * 60)
    print("ПРОВЕРКА РЕШЕНИЯ")
    print("=" * 60)

    print("\n1. Проверка решения A*x = b:")
    max_error = 0
    for i in range(n):
        Ax_i = sum(A[i][j] * x[j] for j in range(n))
        error = abs(Ax_i - b[i])
        max_error = max(max_error, error)
        print(f"  Уравнение {i}: A[{i}]*x = {Ax_i:.8f}, b[{i}] = {b[i]}, ошибка = {error:.2e}")

    print(f"Максимальная ошибка: {max_error:.2e}")
    if max_error < 1e-10:
        print("✅ Решение корректное")
    else:
        print("⚠️ Возможная ошибка в решении")

    print("\n2. Проверка обратной матрицы A*A⁻¹ = I:")
    max_inv_error = 0
    for i in range(n):
        for j in range(n):
            element = sum(A[i][k] * inv[k][j] for k in range(n))
            if i == j:
                error = abs(element - 1)
            else:
                error = abs(element)
            max_inv_error = max(max_inv_error, error)
    
    print(f"Максимальная ошибка в A*A⁻¹: {max_inv_error:.2e}")
    if max_inv_error < 1e-10:
        print("✅ Обратная матрица корректна")
    else:
        print("⚠️ Возможная ошибка в обратной матрице")

    print(f"\n3. Определитель: {det:f}")

    print("\n" + "=" * 60)
    print("РЕЗУЛЬТАТЫ")
    print("=" * 60)

    return x, det, inv


A = [
    [8, 7, 0, 0, 0, 0, 0, 0],
    [7, 13, -4, 0, 0, 0, 0, 0],
    [0, 2, 11, 8, 0, 0, 0, 0],
    [0, 0, -4, 11, -7, 0, 0, 0],
    [0, 0, 0, 6, -13, 5, 0, 0],
    [0, 0, 0, 0, 5, 11, 4, 0],
    [0, 0, 0, 0, 0, 3, 9, 5],
    [0, 0, 0, 0, 0, 0, 4, 5],
    ]
b = [16.0, 26.0, 31.0, 51.0, -53.0, 33.0, 54.0, 35.0]

x, detA, invA = gauss_full(A, b)

print("\nРешение системы:")
print([val for val in x])

print("\nОпределитель:")
print(detA)

print("\nОбратная матрица:")
for row in invA:
    print([val for val in row])

print("\nПроверка A * A⁻¹:")
I_check = [[0.0] * len(A) for _ in range(len(A))]
for i in range(len(A)):
    for j in range(len(A)):
        for k in range(len(A)):
            I_check[i][j] += A[i][k] * invA[k][j]

print("Должна быть единичная матрица:")
for row in I_check:
    print([f"{x:f}" for x in row])