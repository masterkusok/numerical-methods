# i  = |    0   |.  1.   |.  2.   |.   3.  |.   4.  |.   5.  |.   6.  |.   7.   | . .  8  |
# xi = | -2,27  | ―1,83  | -1,39. | -0,95. | -0,51. | -0,07. | 0,37.  | 0,81.   | 1,25    |
# yi = | 3,2383 | 4,5458 | 3,1494 | 2,8715 | 0,8152 | 1,1248 | 0,4714 | ―1,0643 | ―0,1375 |  

import matplotlib.pyplot as plt

xi_full = [-2.27, -1.83, -1.39, -0.95, -0.51, -0.07, 0.37, 0.81, 1.25]
yi_full = [3.2383, 4.5458, 3.1494, 2.8715, 0.8152, 1.1248, 0.4714, -1.0643, -0.1375]
x_star = 0.176


def divided_differences(x_nodes, y_nodes):
    n = len(x_nodes)
    diffs = [y_nodes[:]] 

    for k in range(1, n):
        current_diffs = []
        for i in range(n - k):
            numerator = diffs[k - 1][i + 1] - diffs[k - 1][i]
            denominator = x_nodes[i + k] - x_nodes[i]
            current_diffs.append(numerator / denominator)
        diffs.append(current_diffs)

    coefficients = [d[0] for d in diffs]
    return coefficients

def newton_poly_val(x_val, x_nodes, coefficients):
    n = len(x_nodes)
    result = coefficients[n - 1]
    
    for i in range(n - 2, -1, -1):
        result = result * (x_val - x_nodes[i]) + coefficients[i]
        
    return result

def remainder_estimate(x_val, x_nodes_poly, diffs_err_coeffs):
    product_term = 1.0
    for x_i in x_nodes_poly:
        product_term *= (x_val - x_i)

    divided_difference = diffs_err_coeffs[-1]
    
    return product_term * divided_difference

candidates_p2 = [([4, 5, 6], [4, 5, 6, 7]), ([5, 6, 7], [4, 5, 6, 7])]
candidates_p3 = [([3, 4, 5, 6], [3, 4, 5, 6, 7]), ([4, 5, 6, 7], [4, 5, 6, 7, 8]), ([5, 6, 7, 8], [5, 6, 7, 8])]

print("\nПеребор наборов для многочлена степени 2:")
best_error_p2 = float('inf')
idx_p2 = None
results_p2 = []
for idx, idx_err in candidates_p2:
    x_poly = [xi_full[i] for i in idx]
    y_poly = [yi_full[i] for i in idx]
    x_err = [xi_full[i] for i in idx_err]
    y_err = [yi_full[i] for i in idx_err]
    coeffs_err = divided_differences(x_err, y_err)
    coeffs = divided_differences(x_poly, y_poly)
    error = abs(remainder_estimate(x_star, x_poly, coeffs_err))
    results_p2.append((idx, x_poly, y_poly, coeffs))
    print(f"  Узлы {idx}: погрешность = {error:.6e}")
    if error < best_error_p2:
        best_error_p2 = error
        idx_p2 = idx
print(f"Выбран набор {idx_p2} с погрешностью {best_error_p2:.6e}\n")

x_p2 = [xi_full[i] for i in idx_p2]
y_p2 = [yi_full[i] for i in idx_p2]
coeffs_p2 = divided_differences(x_p2, y_p2)
y_star_p2 = newton_poly_val(x_star, x_p2, coeffs_p2)
r2_estimate = best_error_p2

print("Перебор наборов для многочлена степени 3:")
best_error_p3 = float('inf')
idx_p3 = None
results_p3 = []
for idx, idx_err in candidates_p3:
    x_poly = [xi_full[i] for i in idx]
    y_poly = [yi_full[i] for i in idx]
    x_err = [xi_full[i] for i in idx_err]
    y_err = [yi_full[i] for i in idx_err]
    coeffs_err = divided_differences(x_err, y_err)
    coeffs = divided_differences(x_poly, y_poly)
    error = abs(remainder_estimate(x_star, x_poly, coeffs_err))
    results_p3.append((idx, x_poly, y_poly, coeffs))
    print(f"  Узлы {idx}: погрешность = {error:.6e}")
    if error < best_error_p3:
        best_error_p3 = error
        idx_p3 = idx
print(f"Выбран набор {idx_p3} с погрешностью {best_error_p3:.6e}\n")

x_p3 = [xi_full[i] for i in idx_p3]
y_p3 = [yi_full[i] for i in idx_p3]
coeffs_p3 = divided_differences(x_p3, y_p3)
y_star_p3 = newton_poly_val(x_star, x_p3, coeffs_p3)
r3_estimate = best_error_p3

check_node_x_p2 = x_p2[1]
check_node_y_real_p2 = yi_full[idx_p2[1]]
check_node_y_calc_p2 = newton_poly_val(check_node_x_p2, x_p2, coeffs_p2)

check_node_x_p3 = x_p3[2]
check_node_y_real_p3 = yi_full[idx_p3[2]]
check_node_y_calc_p3 = newton_poly_val(check_node_x_p3, x_p3, coeffs_p3)

print(f"Точка интерполяции x*: {x_star}")

print("Многочлен степени 2:")
print(f"Узлы: {x_p2}")
print(f"Значение: {y_star_p2:.6f}")
print(f"Оценка погрешности: {abs(r2_estimate):.6e}")
print(f"Проверка P₂(x) в узле x={check_node_x_p2:.2f}:")
print(f"  Табличное значение y: {check_node_y_real_p2:.6f}")
print(f"  Вычисленное по P₂:    {check_node_y_calc_p2:.6f}")
print("")

print("Многочлен степени 3:")
print(f"Узлы: {x_p3}")
print(f"Значение: {y_star_p3:.6f}")
print(f"Оценка погрешности: {abs(r3_estimate):.6e}")
print(f"Проверка P₃(x) в узле x={check_node_x_p3:.2f}:")
print(f"  Табличное значение y: {check_node_y_real_p3:.6f}")
print(f"  Вычисленное по P₃:    {check_node_y_calc_p3:.6f}")


def generate_x_range(start, end, num_points):
    step = (end - start) / (num_points - 1)
    return [start + i * step for i in range(num_points)]

x_min = min(xi_full) - 0.2
x_max = max(xi_full) + 0.2
x_plot = generate_x_range(x_min, x_max, 100)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(f'Интерполяция многочленом Ньютона в точке x* = {x_star}', fontsize=16)

ax1.plot(xi_full, yi_full, 'o', color='lightgray', label='Исходные узлы')
colors_p2 = ['b-', 'c-']
for i, (idx, x_poly, y_poly, coeffs) in enumerate(results_p2):
    y_plot = [newton_poly_val(x, x_poly, coeffs) for x in x_plot]
    linewidth = 3 if idx == idx_p2 else 1
    label = f'P₂ узлы {idx}' + (' (лучший)' if idx == idx_p2 else '')
    ax1.plot(x_plot, y_plot, colors_p2[i], linewidth=linewidth, label=label)
ax1.plot(x_p2, y_p2, 'ro', label='Узлы лучшего P₂')
ax1.plot(x_star, y_star_p2, 'g*', markersize=15, label=f'P₂({x_star})={y_star_p2:.4f}')

ax1.set_title(f'Многочлен Ньютона 2-й степени (3 узла)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.legend()

ax2.plot(xi_full, yi_full, 'o', color='lightgray', label='Исходные узлы')
colors_p3 = ['m-', 'r-', 'y-']
for i, (idx, x_poly, y_poly, coeffs) in enumerate(results_p3):
    y_plot = [newton_poly_val(x, x_poly, coeffs) for x in x_plot]
    linewidth = 3 if idx == idx_p3 else 1
    label = f'P₃ узлы {idx}' + (' (лучший)' if idx == idx_p3 else '')
    ax2.plot(x_plot, y_plot, colors_p3[i], linewidth=linewidth, label=label)
ax2.plot(x_p3, y_p3, 'ro', label='Узлы лучшего P₃')
ax2.plot(x_star, y_star_p3, 'g*', markersize=15, label=f'P₃({x_star})={y_star_p3:.4f}')

ax2.set_title(f'Многочлен Ньютона 3-й степени (4 узла)')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.grid(True)
ax2.legend()

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()