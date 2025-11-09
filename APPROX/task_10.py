# i  = |    0   |.  1.   |.  2.   |.   3.  |.   4.  |.   5.  |.   6.  |.   7.   | . .  8  |
# xi = | -2,27  | ―1,83  | -1,39. | -0,95. | -0,51. | -0,07. | 0,37.  | 0,81.   | 1,25    |
# yi = | 3,2383 | 4,5458 | 3,1494 | 2,8715 | 0,8152 | 1,1248 | 0,4714 | ―1,0643 | ―0,1375 |  

x = [-2.27, -1.83, -1.39, -0.95, -0.51, -0.07, 0.37, 0.81, 1.25]
y = [3.2383, 4.5458, 3.1494, 2.8715, 0.8152, 1.1248, 0.4714, -1.0643, -0.1375]

def divided_diff(x_data, y_data):
    n = len(x_data)
    table = [[0] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = y_data[i]
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x_data[i + j] - x_data[i])
    return [table[0][i] for i in range(n)]

def newton(x_data, y_data, x_val):
    coef = divided_diff(x_data, y_data)
    result = coef[0]
    prod = 1
    for i in range(1, len(coef)):
        prod *= (x_val - x_data[i - 1])
        result += coef[i] * prod
    return result

def error_estimate(x_data, x_val, next_coef):
    prod = 1
    for xi in x_data:
        prod *= (x_val - xi)
    return abs(next_coef * prod)

x_star = 0.176

distances = [(abs(xi - x_star), i) for i, xi in enumerate(x)]
distances.sort()

indices_2 = sorted([distances[i][1] for i in range(3)])
x2 = [x[i] for i in indices_2]
y2 = [y[i] for i in indices_2]

indices_3 = sorted([distances[i][1] for i in range(4)])
x3 = [x[i] for i in indices_3]
y3 = [y[i] for i in indices_3]

indices_4 = sorted([distances[i][1] for i in range(5)])
x4 = [x[i] for i in indices_4]
y4 = [y[i] for i in indices_4]

print("Интерполяционный многочлен Ньютона 2-й степени")
print(f"Узлы: x = {x2}")
print(f"      y = {y2}")
coef2 = divided_diff(x2, y2)
print(f"Разделенные разности: {[f'{c:.4f}' for c in coef2]}")
check_val = newton(x2, y2, x2[1])
print(f"Проверка в узле x = {x2[1]}: N(x) = {check_val:.4f}, y = {y2[1]:.4f}")
result_2 = newton(x2, y2, x_star)
coef_next = divided_diff(x3, y3)[3]
error_2 = error_estimate(x2, x_star, coef_next)
print(f"Значение в x* = {x_star}: N(x*) = {result_2:.4f}")
print(f"Оценка погрешности: {error_2:.6f}\n")

print("Интерполяционный многочлен Ньютона 3-й степени")
print(f"Узлы: x = {x3}")
print(f"      y = {y3}")
coef3 = divided_diff(x3, y3)
print(f"Разделенные разности: {[f'{c:.4f}' for c in coef3]}")
check_val = newton(x3, y3, x3[1])
print(f"Проверка в узле x = {x3[1]}: N(x) = {check_val:.4f}, y = {y3[1]:.4f}")
result_3 = newton(x3, y3, x_star)
coef_next = divided_diff(x4, y4)[4]
error_3 = error_estimate(x3, x_star, coef_next)
print(f"Значение в x* = {x_star}: N(x*) = {result_3:.4f}")
print(f"Оценка погрешности: {error_3:.6f}")
