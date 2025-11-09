import copy

def ul_solve_in_place(A_in, b_in, eps=1e-12):
    """
    UL-разложение на месте
    """
    A = copy.deepcopy(A_in)
    b = copy.deepcopy(b_in)
    n = len(A)
    
    for i in range(n):
        A[i][i] = A_in[i][i]

    print("=== Этап UL-разложения) ===")
    print("Исходная матрица A:")
    for row in A:
        print([v for v in row])
    
    for k in range(n-1, -1, -1):
        print(f"\n--- Шаг k={k} ---")
        
        for i in range(0, k+1):
            s = 0.0
            for m in range(k+1, n):
                s += A[i][m] * A[m][k]

            A[i][k] = A_in[i][k] - s
        
        for j in range(0, k):
            s = 0.0
            for m in range(k+1, n):
                s += A[k][m] * A[m][j]

            A[k][j] = (A_in[k][j] - s) / A[k][k]
        
        
        print(f"После шага k={k}:")
        print("Матрица [U|L] (верхний треугольник + диагональ = U, нижний = L):")
        for i in range(n):
            row_str = []
            for j in range(n):
                if i <= j:
                    row_str.append(f"U[{i},{j}]={A[i][j]:.3f}")
                else:
                    row_str.append(f"L[{i},{j}]={A[i][j]:.3f}")
            print(row_str)

    det = 1.0
    for i in range(n):
        det *= A[i][i]

    print(f"\nОпределитель: det = {det:.6f}")

    
    y = [0.0] * n
    for i in range(n-1, -1, -1):
        s = 0.0
        for j in range(i+1, n):
            s += A[i][j] * y[j]
        y[i] = (b[i] - s) / A[i][i]
    
    print("\nВектор y (после U * y = b):", [v for v in y])

    x = [0.0] * n
    for i in range(n):
        s = 0.0
        for j in range(i):
            s += A[i][j] * x[j]
        x[i] = y[i] - s
    
    print("Вектор x (решение):", [v for v in x])

    Ax = [sum(A_in[i][j] * x[j] for j in range(n)) for i in range(n)]
    residual = [Ax[i] - b_in[i] for i in range(n)]
    max_error = max(abs(r) for r in residual)
    
    print("\nПроверка решения:")
    print("A * x =", [v for v in Ax])
    print("b =", b_in)
    print("Невязка (A*x - b) =", [v for v in residual])
    print(f"Максимальная ошибка: {max_error:.2e}")
    
    if max_error < eps:
        print("Решение корректное (невязка в пределах точности)")
    else:
        print("Решение может быть неточным!")

    print("\n=== Вычисление обратной матрицы ===")
    inv = [[0.0] * n for _ in range(n)]
    
    for col in range(n):
        e = [0.0] * n
        e[col] = 1.0
        
        ycol = [0.0] * n
        for i in range(n-1, -1, -1):
            s = 0.0
            for j in range(i+1, n):
                s += A[i][j] * ycol[j]
            ycol[i] = (e[i] - s) / A[i][i]
        
        xcol = [0.0] * n
        for i in range(n):
            s = 0.0
            for j in range(i):
                s += A[i][j] * xcol[j]  # L[i][j]
            xcol[i] = ycol[i] - s  # L[i][i] = 1
        
        for i in range(n):
            inv[i][col] = xcol[i]
    
    # Округлённые результаты
    x_rounded = [v for v in x]
    det_rounded = det
    inv_rounded = [[v for v in row] for row in inv]
    
    # Восстанавливаем отдельные L и U для вывода
    L_out = [[0.0] * n for _ in range(n)]
    U_out = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            if i > j:
                L_out[i][j] = A[i][j]
                U_out[i][j] = 0.0
            else:
                L_out[i][j] = 1.0 if i == j else 0.0
                U_out[i][j] = A[i][j]

    return {
        "x": x, "det": det, "inv": inv,
        "x_rounded": x_rounded, "det_rounded": det_rounded, "inv_rounded": inv_rounded,
        "L": L_out, "U": U_out, "UL_combined": A
    }

def print_combined_matrix(A, round_digits=3):
    """
    Красиво печатает объединённую матрицу UL
    """
    n = len(A)
    print("\nОбъединённая матрица UL:")
    print("(Верхний треугольник + диагональ = U, нижний треугольник = L)")
    print("-" * 60)
    
    for i in range(n):
        row_str = ""
        for j in range(n):
            if i < j:  # Верхний треугольник (U)
                row_str += f"U{A[i][j]:>8} "
            elif i == j:  # Диагональ (U)
                row_str += f"D{A[i][j]:>8} "
            else:  # Нижний треугольник (L)
                row_str += f"L{A[i][j]:>8} "
        print(row_str)

if __name__ == "__main__":
    A = [
        [4, -6, -3, 9, 4],
        [-5, 3, 1, 4, 2],
        [-2, 1, 4, 2, -5],
        [7, 5, 0, 1, 3],
        [8, -3, 4, 1, -2]
    ]
    b = [31, 62, 47, 77, 52]

    print("Исходная система:")
    print("A =")
    for row in A:
        print([f"{x}" for x in row])
    print("b =", b)
    
    res = ul_solve_in_place(A, b)
    
    print("\n" + "="*50)
    print("ИТОГИ")
    print("="*50)
    
    print_combined_matrix(res["UL_combined"])
    
    print(f"\nРешение x: {res['x_rounded']}")
    print(f"Определитель: {res['det_rounded']}")
    
    print("\nОбратная матрица A⁻¹:")
    for row in res["inv_rounded"]:
        print([f"{x:f}" for x in row])
    
    # Проверка A * A⁻¹ = I
    print("\nПроверка A * A⁻¹:")
    I_check = [[0.0] * len(A) for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A)):
            for k in range(len(A)):
                I_check[i][j] += A[i][k] * res["inv_rounded"][k][j]
    
    print("Должна быть единичная матрица:")
    for row in I_check:
        print([f"{x:f}" for x in row])