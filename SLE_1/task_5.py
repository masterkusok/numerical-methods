def print_matrix(matrix):
    for row in matrix:
        for element in row:
            print(f"{element:8.2f}", end=" ")
        print()


def print_vector(x, iteration):
    print(f"\nИтерация {iteration}:")
    for i in range(len(x)):
        print(f"x{i + 1} = {x[i]:.6f}")


def check_diagonal_dominance(matrix):
    n = len(matrix)
    is_dominant = True

    print("\nПроверка диагонального преобладания:")
    print("-" * 60)

    for i in range(n):
        diagonal = abs(matrix[i][i])
        sum_others = 0
        for j in range(n):
            if i != j:
                sum_others = sum_others + abs(matrix[i][j])

        print(f"Строка {i + 1}: |a{i + 1}{i + 1}| = {diagonal:.2f}, сумма остальных = {sum_others:.2f}", end="")

        if diagonal >= sum_others:
            print(" ✓")
        else:
            print(" ✗")
            is_dominant = False

    return is_dominant


def rearrange_system(matrix, b):
    n = len(matrix)
    used_rows = []
    new_matrix = []
    new_b = []

    for col in range(n):
        best_row = -1
        best_ratio = 0

        for row in range(n):
            if row in used_rows:
                continue

            diagonal = abs(matrix[row][col])
            sum_others = 0
            for j in range(n):
                if j != col:
                    sum_others = sum_others + abs(matrix[row][j])

            if sum_others > 0:
                ratio = diagonal / sum_others
            else:
                ratio = diagonal

            if ratio > best_ratio:
                best_ratio = ratio
                best_row = row

        if best_row != -1:
            used_rows.append(best_row)
            new_matrix.append(matrix[best_row][:])
            new_b.append(b[best_row])

    return new_matrix, new_b


def seidel_method(matrix, b, epsilon, max_iterations):
    n = len(matrix)
    x = []
    for i in range(n):
        x.append(0.0)

    x_new = []
    for i in range(n):
        x_new.append(0.0)

    print("\nНачальное приближение:")
    for i in range(n):
        print(f"x{i + 1} = {x[i]:.6f}")

    for iteration in range(1, max_iterations + 1):
        for i in range(n):
            sum_val = b[i]

            for j in range(n):
                if j < i:
                    sum_val = sum_val - matrix[i][j] * x_new[j]
                elif j > i:
                    sum_val = sum_val - matrix[i][j] * x[j]

            x_new[i] = sum_val / matrix[i][i]

        max_diff = 0
        for i in range(n):
            diff = abs(x_new[i] - x[i])
            if diff > max_diff:
                max_diff = diff

        for i in range(n):
            x[i] = x_new[i]

        if iteration <= 10 or iteration % 10 == 0:
            print_vector(x, iteration)
            print(f"Максимальное отклонение: {max_diff:.8f}")

        if max_diff < epsilon:
            print_vector(x, iteration)
            print(f"Максимальное отклонение: {max_diff:.8f}")
            print(f"\nСходимость достигнута за {iteration} итераций")
            return x, iteration

    print(f"\nДостигнуто максимальное число итераций: {max_iterations}")
    return x, max_iterations


def check_solution(matrix, b, x):
    n = len(x)
    print("\n" + "-" * 60)
    print("ПРОВЕРКА РЕШЕНИЯ:")
    print("-" * 60)

    for i in range(n):
        summ = 0
        for j in range(n):
            summ = summ + matrix[i][j] * x[j]

        error = abs(summ - b[i])
        print(f"Уравнение {i + 1}: {summ:.6f} ≈ {b[i]:.2f} (погрешность: {error:.8f})")

matrix = [
    [5, 15, -4, -6],
    [6, -7, 3, 18],
    [7, -2, -14, -3],
    [19, -4, 7, -5]
]

b = [72, 53, -50, -20]

print("\nИсходная система уравнений:")
print("-" * 60)
print("\nМатрица коэффициентов A:")
print_matrix(matrix)
print("\nВектор правых частей b:")
for i in range(len(b)):
    print(f"b{i + 1} = {b[i]:.2f}")

is_dominant = check_diagonal_dominance(matrix)

if not is_dominant:
    print("\nДиагональное преобладание не выполнено.")
    print("Выполняем перестановку уравнений...")

    matrix, b = rearrange_system(matrix, b)

    print("\nПреобразованная система:")
    print("-" * 60)
    print("\nМатрица коэффициентов A:")
    print_matrix(matrix)
    print("\nВектор правых частей b:")
    for i in range(len(b)):
        print(f"b{i + 1} = {b[i]:.2f}")

    is_dominant = check_diagonal_dominance(matrix)

    if not is_dominant:
        print("\nВнимание! Даже после перестановки диагональное преобладание слабое.")
        print("Метод может сходиться медленно.")

epsilon = 0.0001
max_iterations = 1000

print("\n" + "-" * 60)
print(f"Параметры метода:")
print(f"Точность (epsilon): {epsilon}")
print(f"Максимальное число итераций: {max_iterations}")
print("-" * 60)

solution, iterations = seidel_method(matrix, b, epsilon, max_iterations)

print("\n" + "-" * 60)
print("ИТОГОВОЕ РЕШЕНИЕ:")
print("-" * 60)
for i in range(len(solution)):
    print(f"x{i + 1} = {solution[i]:.8f}")

print(f"\nКоличество итераций: {iterations}")

check_solution(matrix, b, solution)

print("\n" + "-" * 60)
print("КОНЕЦ РАБОТЫ")
print("-" * 60)