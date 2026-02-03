import matplotlib.pyplot as plt
import math

xi = [-2.27, -1.83, -1.39, -0.95, -0.51, -0.07, 0.37, 0.81, 1.25]
yi = [3.2383, 4.5458, 3.1494, 2.8715, 0.8152, 1.1248, 0.4714, -1.0643, -0.1375]
x_star = 0.176

def lagrange_interpolation(nodes_x, nodes_y, x_val):
    n = len(nodes_x)
    Ln = 0.0
    for i in range(n):
        l_i = 1.0
        for j in range(n):
            if i != j:
                l_i *= (x_val - nodes_x[j]) / (nodes_x[i] - nodes_x[j])
        Ln += nodes_y[i] * l_i
    return Ln

indices_L2 = [5, 6, 7]
indices_L2_2 = [4, 5, 6]
all_indices_L2 = [indices_L2, indices_L2_2]

indices_L3 = [4, 5, 6, 7]
indices_L3_2 = [3, 4, 5, 6]
indices_L3_3 = [5, 6, 7, 8]
all_indices_L3 = [indices_L3, indices_L3_2, indices_L3_3]

nodes_x_2 = [xi[i] for i in indices_L2]
nodes_y_2 = [yi[i] for i in indices_L2]
nodes_x_3 = [xi[i] for i in indices_L3]
nodes_y_3 = [yi[i] for i in indices_L3]

val_L2 = lagrange_interpolation(nodes_x_2, nodes_y_2, x_star)
val_L3 = lagrange_interpolation(nodes_x_3, nodes_y_3, x_star)

def generate_x_range(start, end, num_points):
    step = (end - start) / (num_points - 1)
    return [start + i * step for i in range(num_points)]

x_plot = generate_x_range(min(xi)-0.2, max(xi)+0.2, 50)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.plot(xi, yi, 'o', color='lightgray', label='Остальные узлы')
colors_L2 = ['b-', 'c-']
for i, idx in enumerate(all_indices_L2):
    nodes_x = [xi[j] for j in idx]
    nodes_y = [yi[j] for j in idx]
    y_plot = [lagrange_interpolation(nodes_x, nodes_y, x) for x in x_plot]
    linewidth = 3 if idx == indices_L2 else 1
    label = f'L2 узлы {idx}' + (' (основной)' if idx == indices_L2 else '')
    ax1.plot(x_plot, y_plot, colors_L2[i], linewidth=linewidth, label=label)
ax1.plot(nodes_x_2, nodes_y_2, 'ro', label='Узлы основного L2')
ax1.plot(x_star, val_L2, 'g*', markersize=15, label=f'x*={x_star}')

ax1.set_title(f'Лагранж 2-й степени (3 узла)\ny(x*) = {val_L2:.4f}')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.legend()

ax2.plot(xi, yi, 'o', color='lightgray', label='Остальные узлы')
colors_L3 = ['m-', 'r-', 'y-']
for i, idx in enumerate(all_indices_L3):
    nodes_x = [xi[j] for j in idx]
    nodes_y = [yi[j] for j in idx]
    y_plot = [lagrange_interpolation(nodes_x, nodes_y, x) for x in x_plot]
    linewidth = 3 if idx == indices_L3 else 1
    label = f'L3 узлы {idx}' + (' (основной)' if idx == indices_L3 else '')
    ax2.plot(x_plot, y_plot, colors_L3[i], linewidth=linewidth, label=label)
ax2.plot(nodes_x_3, nodes_y_3, 'ro', label='Узлы основного L3')
ax2.plot(x_star, val_L3, 'g*', markersize=15, label=f'x*={x_star}')

ax2.set_title(f'Лагранж 3-й степени (4 узла)\ny(x*) = {val_L3:.4f}')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()

print(f"Значение функции в точке {x_star} (Степень 2): {val_L2:.5f}")
print(f"Значение функции в точке x={nodes_x_2[0]} (Степень 2): {lagrange_interpolation(nodes_x_2, nodes_y_2, nodes_x_2[0])}")
print(f"Значение функции в точке {x_star} (Степень 3): {val_L3:.5f}")
print(f"Значение функции в точке x={nodes_x_3[0]} (Степень 3): {lagrange_interpolation(nodes_x_3, nodes_y_3, nodes_x_3[0])}")
