import math
import matplotlib.pyplot as plt

def F(x):
    return math.sin(x) * math.cosh(2*x) / (x*x+1)

a = 0
b = 3.2

def midpoint_rule(f, a, b, n):
    h = (b - a) / n
    s = sum(f(a + (i + 0.5) * h) for i in range(n))
    return h * s

def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b)) + sum(f(a + i * h) for i in range(1, n))
    return h * s

def simpson_rule(f, a, b, n):
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4 * sum(f(a + i * h) for i in range(1, n, 2))
    s += 2 * sum(f(a + i * h) for i in range(2, n, 2))
    return h * s / 3

def dF(x):
    s = math.sin(x)
    c = math.cos(x)
    sh = math.sinh(2*x)
    ch = math.cosh(2*x)
    x2_1 = x**2 + 1
    
    numerator = (c * ch + 2 * s * sh) * x2_1 - (2 * x * s * ch)
    denominator = x2_1**2
    
    return numerator / denominator

def euler_rule(f, a, b, n):
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b)) + sum(f(a + i * h) for i in range(1, n))
    res = h * s
    correction = (h**2 / 12) * (dF(a) - dF(b))
    return res + correction

def runge_romberg(I_h, I_h2, p):
    return I_h2 + (I_h2 - I_h) / (2**p - 1)

def exact_integral(f, a, b, n=10000):
    return trapezoidal_rule(f, a, b, n)

n = 8
n2 = n * 2

I_exact = exact_integral(F, a, b)

print(f"Точное значение интеграла: {I_exact:.10f}\n")

methods = [
    ("Средних прямоугольников", midpoint_rule, 2),
    ("Трапеций", trapezoidal_rule, 2),
    ("Симпсона", simpson_rule, 4),
    ("Эйлера", euler_rule, 4)
]

results = []

for name, method, p in methods:
    I_h = method(F, a, b, n)
    I_h2 = method(F, a, b, n2)
    I_rr = runge_romberg(I_h, I_h2, p)
    
    err_h = abs(I_h - I_exact)
    err_h2 = abs(I_h2 - I_exact)
    err_rr = abs(I_rr - I_exact)
    
    results.append((name, I_h, I_h2, I_rr, err_h, err_h2, err_rr))
    
    print(f"Метод {name}:")
    print(f"  n = {n}:  I = {I_h:.10f}, погрешность = {err_h:.2e}")
    print(f"  n = {n2}: I = {I_h2:.10f}, погрешность = {err_h2:.2e}")
    print(f"  Рунге-Ромберг: I = {I_rr:.10f}, погрешность = {err_rr:.2e}")
    print()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

method_names = [r[0] for r in results]
errors_h = [r[4] for r in results]
errors_h2 = [r[5] for r in results]
errors_rr = [r[6] for r in results]

x_pos = range(len(method_names))
width = 0.25

ax1.bar([x - width for x in x_pos], errors_h, width, label=f'n = {n}', alpha=0.8)
ax1.bar(x_pos, errors_h2, width, label=f'n = {n2}', alpha=0.8)
ax1.bar([x + width for x in x_pos], errors_rr, width, label='Рунге-Ромберг', alpha=0.8)
ax1.set_yscale('log')
ax1.set_ylabel('Погрешность')
ax1.set_title('Сравнение погрешностей методов')
ax1.set_xticks(x_pos)
ax1.set_xticklabels(method_names, rotation=15, ha='right')
ax1.legend()
ax1.grid(True, alpha=0.3)

x_vals = [a + i * (b - a) / 100 for i in range(101)]
y_vals = [F(x) for x in x_vals]
ax2.plot(x_vals, y_vals, 'b-', linewidth=2)
ax2.fill_between(x_vals, 0, y_vals, alpha=0.3)
ax2.set_xlabel('x')
ax2.set_ylabel('F(x)')
ax2.set_title('Подынтегральная функция')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

fig2, axes = plt.subplots(2, 2, figsize=(14, 10))

n_vis = 8
h_vis = (b - a) / n_vis
x_points = [a + i * h_vis for i in range(n_vis + 1)]

ax = axes[0, 0]
x_fine = [a + i * (b - a) / 200 for i in range(201)]
y_fine = [F(x) for x in x_fine]
ax.plot(x_fine, y_fine, 'r-', linewidth=2, label='F(x)')
for i in range(n_vis):
    x_left = a + i * h_vis
    x_right = a + (i + 1) * h_vis
    x_mid = x_left + h_vis / 2
    y_mid = F(x_mid)
    ax.bar(x_mid, y_mid, width=h_vis, alpha=0.3, edgecolor='blue', linewidth=1.5)
    ax.plot([x_left, x_right], [y_mid, y_mid], 'b-', linewidth=1.5)
ax.set_xlabel('x')
ax.set_ylabel('F(x)')
ax.set_title(f'Метод средних прямоугольников\nn={n}: {results[0][1]:.6f}, n={n2}: {results[0][2]:.6f}, РР: {results[0][3]:.6f}', fontsize=10)
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(x_fine, y_fine, 'r-', linewidth=2, label='F(x)')
for i in range(n_vis):
    x_left = x_points[i]
    x_right = x_points[i + 1]
    y_left = F(x_left)
    y_right = F(x_right)
    ax.fill([x_left, x_right, x_right, x_left], [0, 0, y_right, y_left], alpha=0.3, color='green')
    ax.plot([x_left, x_right], [y_left, y_right], 'g-', linewidth=1.5)
ax.set_xlabel('x')
ax.set_ylabel('F(x)')
ax.set_title(f'Метод трапеций\nn={n}: {results[1][1]:.6f}, n={n2}: {results[1][2]:.6f}, РР: {results[1][3]:.6f}', fontsize=10)
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.plot(x_fine, y_fine, 'r-', linewidth=2, label='F(x)')
for i in range(0, n_vis, 2):
    if i + 2 <= n_vis:
        x0, x1, x2 = x_points[i], x_points[i + 1], x_points[i + 2]
        y0, y1, y2 = F(x0), F(x1), F(x2)
        x_parab = [x0 + j * (x2 - x0) / 50 for j in range(51)]
        y_parab = []
        for x in x_parab:
            L0 = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2))
            L1 = ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2))
            L2 = ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1))
            y_parab.append(y0 * L0 + y1 * L1 + y2 * L2)
        ax.fill_between(x_parab, 0, y_parab, alpha=0.3, color='orange')
        ax.plot(x_parab, y_parab, 'orange', linewidth=1.5)
ax.set_xlabel('x')
ax.set_ylabel('F(x)')
ax.set_title(f'Метод Симпсона (параболы)\nn={n}: {results[2][1]:.6f}, n={n2}: {results[2][2]:.6f}, РР: {results[2][3]:.6f}', fontsize=10)
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(x_fine, y_fine, 'r-', linewidth=2, label='F(x)')
for i in range(n_vis):
    x_left = x_points[i]
    x_right = x_points[i + 1]
    y_left = F(x_left)
    y_right = F(x_right)
    ax.fill([x_left, x_right, x_right, x_left], [0, 0, y_right, y_left], alpha=0.3, color='purple')
    ax.plot([x_left, x_right], [y_left, y_right], color='purple', linewidth=1.5)
for x_pt in x_points:
    ax.plot(x_pt, F(x_pt), 'ko', markersize=5)
ax.set_xlabel('x')
ax.set_ylabel('F(x)')
ax.set_title(f'Метод Эйлера\nn={n}: {results[3][1]:.6f}, n={n2}: {results[3][2]:.6f}, РР: {results[3][3]:.6f}', fontsize=10)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()