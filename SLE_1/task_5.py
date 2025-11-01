def print_matrix(matrix):
    for row in matrix:
        for element in row:
            print(f"{element:8.4f}", end=" ")
        print()


def print_vector(x, iteration):
    print(f"\nИтерация {iteration}:")
    for i in range(len(x)):
        print(f"x{i + 1} = {x[i]:.8f}")


def transform_to_iteration_form(matrix, b):
    n = len(matrix)
    alpha = []
    beta = []

    for i in range(n):
        row_alpha = []
        for j in range(n):
            if i == j:
                row_alpha.append(0.0)
            else:
                row_alpha.append(-matrix[i][j] / matrix[i][i])
        alpha.append(row_alpha)
        beta.append(b[i] / matrix[i][i])

    return alpha, beta


def check_convergence(alpha):
    n = len(alpha)
    max_row_sum = 0.0
    max_col_sum = 0.0

    print("\nПроверка условия сходимости:")
    for i in range(n):
        row_sum = sum(abs(alpha[i][j]) for j in range(n))
        if row_sum > max_row_sum:
            max_row_sum = row_sum

    print("\n▶ Проверка по столбцам (1-норма):")
    for j in range(n):
        col_sum = sum(abs(alpha[i][j]) for i in range(n))
        if col_sum > max_col_sum:
            max_col_sum = col_sum

    if max_row_sum < 1 or max_col_sum < 1:
        print("\nУсловие сходимости выполнено (одна из норм < 1)")
        return True
    else:
        print("\Условие сходимости не выполнено (все нормы >= 1)")
        return False


def zeidel_method(alpha, beta, epsilon, max_iterations):
    n = len(alpha)
    x = beta[:]  # начальное приближение

    print("\nНачальное приближение x(0) = beta:")
    for i in range(n):
        print(f"x{i + 1} = {x[i]:.8f}")

    for iteration in range(1, max_iterations + 1):
        x_old = x[:]
        max_diff = 0.0

        for i in range(n):
            sum_val = beta[i]
            for j in range(n):
                if j != i:
                    if j < i:
                        sum_val += alpha[i][j] * x[j]
                    else:
                        sum_val += alpha[i][j] * x_old[j] 
            x[i] = sum_val

            diff = abs(x[i] - x_old[i])
            if diff > max_diff:
                max_diff = diff

        if iteration <= 10 or iteration % 10 == 0:
            print_vector(x, iteration)
            print(f"Максимальное отклонение: {max_diff:.10f}")

        if max_diff < epsilon:
            print_vector(x, iteration)
            print(f"Максимальное отклонение: {max_diff:.10f}")
            print(f"\nСходимость достигнута за {iteration} итераций ✓")
            return x, iteration

    print(f"\nДостигнуто максимальное число итераций: {max_iterations}")
    return x, max_iterations



def check_solution(matrix, b, x):
    n = len(x)
    print("ПРОВЕРКА РЕШЕНИЯ:")

    for i in range(n):
        summ = 0
        for j in range(n):
            summ = summ + matrix[i][j] * x[j]

        error = abs(summ - b[i])
        print(f"Уравнение {i + 1}: {summ:.8f} ≈ {b[i]:.2f} (погрешность: {error:.10f})")


def rearrange_for_convergence(matrix, b):
    n = len(matrix)
    used_rows = []
    new_matrix = []
    new_b = []

    for col in range(n):
        best_row = -1
        best_value = 0

        for row in range(n):
            if row in used_rows:
                continue

            if abs(matrix[row][col]) > best_value:
                best_value = abs(matrix[row][col])
                best_row = row

        if best_row != -1:
            used_rows.append(best_row)
            new_matrix.append(matrix[best_row][:])
            new_b.append(b[best_row])

    return new_matrix, new_b

matrix = [
    [5, 15, -4, -6],
    [6, -7, 3, 18],
    [7, -2, -14, -3],
    [19, -4, 7, -5]
]

b = [72, 53, -50, -20]

print("\nИсходная система уравнений:")
print("\nМатрица коэффициентов A:")
print_matrix(matrix)
print("\nВектор правых частей b:")
for i in range(len(b)):
    print(f"b{i + 1} = {b[i]:.2f}")
print("ПРЕОБРАЗОВАНИЕ К ВИДУ x = alpha*x + beta")

print("\nПерестановка уравнений для улучшения сходимости...")
matrix, b = rearrange_for_convergence(matrix, b)

print("\nПреобразованная система:")
print("\nМатрица коэффициентов A:")
print_matrix(matrix)
print("\nВектор правых частей b:")
for i in range(len(b)):
    print(f"b{i + 1} = {b[i]:.2f}")

alpha, beta = transform_to_iteration_form(matrix, b)

print("\nМатрица alpha:")
print_matrix(alpha)
print("\nВектор beta:")
for i in range(len(beta)):
    print(f"beta{i + 1} = {beta[i]:.8f}")

converges = check_convergence(alpha)

if not converges:
    print("\nВнимание! Условие сходимости не выполнено строго.")
    print("Метод может не сойтись или сходиться медленно.")

epsilon = 0.0001
max_iterations = 1000

print("ИТЕРАЦИОННЫЙ ПРОЦЕСС")
print(f"Точность (epsilon): {epsilon}")
print(f"Максимальное число итераций: {max_iterations}")

solution, iterations = zeidel_method(alpha, beta, epsilon, max_iterations)

print("ИТОГОВОЕ РЕШЕНИЕ:")
for i in range(len(solution)):
    print(f"x{i + 1} = {solution[i]:.10f}")

print(f"\nКоличество итераций: {iterations}")

check_solution(matrix, b, solution)
