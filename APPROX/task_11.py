# **Вариант 41** (x* = 3,074)
# 
# | i | x_i   | y_i   |
# |---|-------|-------|
# | 1 | 5.267 | 3.682 |
# | 2 | 3.928 | 4.298 |
# | 3 | 9.702 | 3.554 |
# | 4 | 5.704 | 6.008 |

import matplotlib.pyplot as plt

data = [(5.267, 3.682), (3.928, 4.298), (9.702, 3.554), (5.704, 6.008)]
data.sort()
x = [p[0] for p in data]
y = [p[1] for p in data]
n = len(x) - 1
x_star = 3.074

# Вычисление h_i
h = [x[i+1] - x[i] for i in range(n)]

# Система для нахождения c_i (метод прогонки)
# Естественные граничные условия: c_0 = c_n = 0
a = [0] * (n + 1)
b = [1] + [2 * (h[i-1] + h[i]) for i in range(1, n)] + [1]
c = [0] + [h[i] for i in range(1, n)] + [0]
d = [0] + [3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]) for i in range(1, n)] + [0]

for i in range(1, n):
    a[i] = h[i-1]

# Прогонка
alpha = [0] * (n + 1)
beta = [0] * (n + 1)
for i in range(1, n + 1):
    alpha[i] = c[i] / (b[i] - a[i] * alpha[i-1])
    beta[i] = (d[i] + a[i] * beta[i-1]) / (b[i] - a[i] * alpha[i-1])

coef_c = [0] * (n + 1)
for i in range(n - 1, -1, -1):
    coef_c[i] = alpha[i] * coef_c[i+1] + beta[i]

# Вычисление коэффициентов a, b, d
coef_a = y[:]
coef_b = [(y[i+1] - y[i]) / h[i] - h[i] * (coef_c[i+1] + 2 * coef_c[i]) / 3 for i in range(n)]
coef_d = [(coef_c[i+1] - coef_c[i]) / (3 * h[i]) for i in range(n)]

# Поиск отрезка для x*
idx = 0
for i in range(n):
    if x[i] <= x_star <= x[i+1]:
        idx = i
        break

# Вычисление значения сплайна в x*
dx = x_star - x[idx]
S_x_star = coef_a[idx] + coef_b[idx] * dx + coef_c[idx] * dx**2 + coef_d[idx] * dx**3

print(f"Значение сплайна в точке x* = {x_star}: S(x*) = {S_x_star:.6f}")
print(f"\nКоэффициенты сплайна на отрезке [{x[idx]}, {x[idx+1]}]:")
print(f"a = {coef_a[idx]:.6f}")
print(f"b = {coef_b[idx]:.6f}")
print(f"c = {coef_c[idx]:.6f}")
print(f"d = {coef_d[idx]:.6f}")

# Графики
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
for i in range(n):
    x_seg = [x[i] + j * h[i] / 50 for j in range(51)]
    y_seg = [coef_a[i] + coef_b[i] * (xi - x[i]) + coef_c[i] * (xi - x[i])**2 + coef_d[i] * (xi - x[i])**3 for xi in x_seg]
    plt.plot(x_seg, y_seg, 'b-')
plt.plot(x, y, 'ro', label='Узлы интерполяции', markersize=8)
plt.plot(x_star, S_x_star, 'g*', markersize=15, label=f'x* = {x_star}')
plt.grid(True)
plt.legend()
plt.title('Кубический сплайн')

plt.subplot(1, 2, 2)
for i in range(n):
    x_seg = [x[i] + j * h[i] / 50 for j in range(51)]
    y_seg = [coef_a[i] + coef_b[i] * (xi - x[i]) + coef_c[i] * (xi - x[i])**2 + coef_d[i] * (xi - x[i])**3 for xi in x_seg]
    plt.plot(x_seg, y_seg, 'b-', label=f'Отрезок {i+1}' if i < n else '')
plt.plot(x, y, 'ro', markersize=8)
for i in range(n):
    plt.axvline(x[i], color='gray', linestyle='--', alpha=0.3)
plt.axvline(x[n], color='gray', linestyle='--', alpha=0.3)
plt.grid(True)
plt.title('Сплайн по отрезкам')

plt.tight_layout()
plt.show()