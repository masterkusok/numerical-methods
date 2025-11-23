# i  = |    0   |.  1.   |.  2.   |.   3.  |.   4.  |.   5.  |.   6.  |.   7.   | . .  8  |
# xi = | -2,27  | ―1,83  | -1,39. | -0,95. | -0,51. | -0,07. | 0,37.  | 0,81.   | 1,25    |
# yi = | 3,2383 | 4,5458 | 3,1494 | 2,8715 | 0,8152 | 1,1248 | 0,4714 | ―1,0643 | ―0,1375 |  

import matplotlib.pyplot as plt

x = [-2.27, -1.83, -1.39, -0.95, -0.51, -0.07, 0.37, 0.81, 1.25]
y = [3.2383, 4.5458, 3.1494, 2.8715, 0.8152, 1.1248, 0.4714, -1.0643, -0.1375]

def lagrange(x_data, y_data, x_val):
    result = 0
    n = len(x_data)
    for i in range(n):
        l = y_data[i]
        for j in range(n):
            if i != j:
                l *= (x_val - x_data[j]) / (x_data[i] - x_data[j])
        result += l
    return result

x_star = 0.176

# Находим ближайшие точки к x*
distances = [(abs(xi - x_star), i) for i, xi in enumerate(x)]
distances.sort()

# Для многочлена 2-й степени (3 точки)
indices_2 = sorted([distances[i][1] for i in range(3)])
x2 = [x[i] for i in indices_2]
y2 = [y[i] for i in indices_2]

# Для многочлена 3-й степени (4 точки)
indices_3 = sorted([distances[i][1] for i in range(4)])
x3 = [x[i] for i in indices_3]
y3 = [y[i] for i in indices_3]

print("Интерполяционный многочлен Лагранжа 2-й степени")
print(f"Узлы: x = {x2}")
print(f"      y = {y2}")
check_node = x2[1]
check_val = lagrange(x2, y2, check_node)
print(f"Проверка в узле x = {check_node}: L(x) = {check_val:.4f}, y = {y2[1]:.4f}")
result_2 = lagrange(x2, y2, x_star)
print(f"Значение в x* = {x_star}: L(x*) = {result_2:.4f}\n")

print("Интерполяционный многочлен Лагранжа 3-й степени")
print(f"Узлы: x = {x3}")
print(f"      y = {y3}")
check_node = x3[1]
check_val = lagrange(x3, y3, check_node)
print(f"Проверка в узле x = {check_node}: L(x) = {check_val:.4f}, y = {y3[1]:.4f}")
result_3 = lagrange(x3, y3, x_star)
print(f"Значение в x* = {x_star}: L(x*) = {result_3:.4f}")

# Графики
x_min, x_max = min(x), max(x)
x_plot = [x_min + i * (x_max - x_min) / 199 for i in range(200)]
y_plot_2 = [lagrange(x2, y2, xi) for xi in x_plot]
y_plot_3 = [lagrange(x3, y3, xi) for xi in x_plot]

plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.plot(x_plot, y_plot_2, 'b-', label='Многочлен 2-й степени')
plt.plot(x2, y2, 'ro', label='Узлы интерполяции')
plt.plot(x_star, result_2, 'g*', markersize=12, label=f'x* = {x_star}')
plt.grid(True)
plt.legend()
plt.title('Многочлен Лагранжа 2-й степени')

plt.subplot(1, 3, 2)
plt.plot(x_plot, y_plot_3, 'r-', label='Многочлен 3-й степени')
plt.plot(x3, y3, 'bo', label='Узлы интерполяции')
plt.plot(x_star, result_3, 'g*', markersize=12, label=f'x* = {x_star}')
plt.grid(True)
plt.legend()
plt.title('Многочлен Лагранжа 3-й степени')

plt.subplot(1, 3, 3)
errors_2 = [abs(lagrange(x2, y2, xi) - lagrange(x, y, xi)) for xi in x_plot]
errors_3 = [abs(lagrange(x3, y3, xi) - lagrange(x, y, xi)) for xi in x_plot]
plt.plot(x_plot, errors_2, 'b-', label='Погрешность 2-й степени')
plt.plot(x_plot, errors_3, 'r-', label='Погрешность 3-й степени')
plt.grid(True)
plt.legend()
plt.title('Сравнение точности')
plt.yscale('log')

plt.tight_layout()
plt.show()
