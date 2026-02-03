import matplotlib.pyplot as plt

xi = [2.30, 2.522, 2.818, 3.262, 3.632, 3.928, 4.298, 4.742, 5.334, 5.704, 6.00]
yi = [2.152, 3.751, 3.394, 1.979, 1.416, 1.795, 2.738, 1.593, 0.468, 0.327, 1.854]
x_star = 3.704

def cubic_spline_interpolation(x_nodes, y_nodes, x_eval):
    n = len(x_nodes) - 1
    h = [x_nodes[i+1] - x_nodes[i] for i in range(n)] 
    
    P = [0.0] * n
    Q = [0.0] * n
    
    for i in range(1, n):
        hi_prev = h[i-1]
        hi_curr = h[i]
        
        A = hi_prev
        B = 2 * (hi_prev + hi_curr)
        C = hi_curr
        
        F = 3 * ((y_nodes[i+1] - y_nodes[i])/hi_curr - (y_nodes[i] - y_nodes[i-1])/hi_prev)
        
        denominator = B + A * P[i-1]
        P[i] = -C / denominator
        Q[i] = (F - A * Q[i-1]) / denominator
        
    c = [0.0] * (n + 1)
    c[n] = 0.0
    
    for i in range(n - 1, 0, -1):
        c[i] = P[i] * c[i+1] + Q[i]
        
    splines = []
    
    target_spline = None
    target_idx = -1

    for i in range(n):
        a = y_nodes[i]
        c_curr = c[i]
        b = (y_nodes[i+1] - y_nodes[i]) / h[i] - (h[i] / 3.0) * (c[i+1] + 2 * c_curr)
        d = (c[i+1] - c_curr) / (3.0 * h[i])
        
        splines.append({'a': a, 'b': b, 'c': c_curr, 'd': d, 'x_start': x_nodes[i]})
        
        if x_nodes[i] <= x_eval <= x_nodes[i+1]:
            target_spline = splines[-1]
            target_idx = i + 1
    
    if target_spline:
        dx = x_eval - target_spline['x_start']
        val = target_spline['a'] + target_spline['b'] * dx + \
              target_spline['c'] * (dx**2) + target_spline['d'] * (dx**3)
        return val, target_spline, target_idx, splines
    else:
        return None, None, None, splines

val_star, coeffs_star, idx_star, all_splines = cubic_spline_interpolation(xi, yi, x_star)

print(f"\nРЕЗУЛЬТАТЫ ДЛЯ ТОЧКИ x* = {x_star}")
print("-" * 50)
if coeffs_star:
    print(f"Точка попадает в интервал {idx_star}: [{xi[idx_star-1]}; {xi[idx_star]}]")
    print(f"Коэффициенты сплайна на этом отрезке:")
    print(f"a = {coeffs_star['a']:.6f}")
    print(f"b = {coeffs_star['b']:.6f}")
    print(f"c = {coeffs_star['c']:.6f}")
    print(f"d = {coeffs_star['d']:.6f}")
    print("-" * 50)
    print(f"Значение функции S(x*) = {val_star:.6f}")
    print("-" * 50)

x_plot = []
y_plot = []
y_prime = []
y_double_prime = []

for i, s in enumerate(all_splines):
    steps = 20
    h_local = xi[i+1] - xi[i]
    for k in range(steps):
        xx = xi[i] + k * h_local / steps
        dx = xx - xi[i]
        yy = s['a'] + s['b']*dx + s['c']*(dx**2) + s['d']*(dx**3)
        yy_prime = s['b'] + 2*s['c']*dx + 3*s['d']*(dx**2)
        yy_double_prime = 2*s['c'] + 6*s['d']*dx
        x_plot.append(xx)
        y_plot.append(yy)
        y_prime.append(yy_prime)
        y_double_prime.append(yy_double_prime)

x_plot.append(xi[-1])
y_plot.append(yi[-1])
dx_last = xi[-1] - xi[-2]
y_prime.append(all_splines[-1]['b'] + 2*all_splines[-1]['c']*dx_last + 3*all_splines[-1]['d']*(dx_last**2))
y_double_prime.append(2*all_splines[-1]['c'] + 6*all_splines[-1]['d']*dx_last)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

ax1.plot(x_plot, y_plot, 'b-', label='Кубический сплайн S(x)')
ax1.plot(xi, yi, 'ro', label='Узлы интерполяции')
ax1.plot(x_star, val_star, 'g*', markersize=12, label=f'Точка x*={x_star}')
ax1.grid(True)
ax1.set_title(f'Интерполяция кубическим сплайном (x*={x_star})')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.legend()

ax2.plot(x_plot, y_prime, 'r-', label="Первая производная S'(x)")
ax2.grid(True)
ax2.set_title("Первая производная сплайна")
ax2.set_xlabel('x')
ax2.set_ylabel("S'(x)")
ax2.legend()

ax3.plot(x_plot, y_double_prime, 'g-', label='Вторая производная S"(x)')
ax3.grid(True)
ax3.set_title('Вторая производная сплайна')
ax3.set_xlabel('x')
ax3.set_ylabel('S"(x)')
ax3.legend()

plt.tight_layout()
plt.show()