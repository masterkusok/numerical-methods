import math
import matplotlib.pyplot as plt

x_start, x_end = 2.5, 9.0
h = 0.5
n = int((x_end - x_start) / h) + 1
x = [x_start + i * h for i in range(n)]

def y_exact(x_val):
    return math.sqrt(x_val * (x_val - 2)) + math.sqrt(x_val)

def p(x_val):
    return -x_val * (x_val - 4) / (2 * x_val**2 * (x_val - 2))

def q(x_val):
    return (x_val - 3) / (2 * x_val**2 * (x_val - 2))

def solve_tridiagonal(a, b, c, d):
    n = len(d)
    alpha, beta = [0] * n, [0] * n
    alpha[0] = -c[0] / b[0]
    beta[0] = d[0] / b[0]
    
    for i in range(1, n):
        denom = b[i] + a[i] * alpha[i-1]
        alpha[i] = -c[i] / denom if i < n-1 else 0
        beta[i] = (d[i] - a[i] * beta[i-1]) / denom
    
    y = [0] * n
    y[-1] = beta[-1]
    for i in range(n-2, -1, -1):
        y[i] = alpha[i] * y[i+1] + beta[i]
    
    return y

def solve_first_order():
    a, b, c, d = [0] * n, [0] * n, [0] * n, [0] * n
    
    b[0] = 3 - 1/h
    c[0] = 1/h
    d[0] = 9.755
    
    for i in range(1, n-1):
        a[i] = 1/h**2 - p(x[i])/(2*h)
        b[i] = -2/h**2 + q(x[i])
        c[i] = 1/h**2 + p(x[i])/(2*h)
        d[i] = 0
    
    b[n-1] = 1
    d[n-1] = 10.937
    
    return solve_tridiagonal(a, b, c, d)

def solve_second_order():
    a, b, c, d = [0] * n, [0] * n, [0] * n, [0] * n
    
    b[0] = 3 + 3/(2*h)
    c[0] = -4/(2*h)
    d[0] = 9.755
    d[0] += 1/(2*h) * 0
    
    a[1] = 1/h**2 - p(x[1])/(2*h)
    b[1] = -2/h**2 + q(x[1])
    c[1] = 1/h**2 + p(x[1])/(2*h)
    d[1] = -a[1] * (1/(2*h)) / (3 + 3/(2*h)) * 9.755
    
    b[0] = 3 + 3/(2*h)
    c[0] = -2/h
    d[0] = 9.755
    
    for i in range(1, n-1):
        a[i] = 1/h**2 - p(x[i])/(2*h)
        b[i] = -2/h**2 + q(x[i])
        c[i] = 1/h**2 + p(x[i])/(2*h)
        d[i] = 0
    
    b[n-1] = 1
    d[n-1] = 10.937
    
    return solve_tridiagonal(a, b, c, d)

y1 = solve_first_order()
y2 = solve_second_order()
y_an = [y_exact(xi) for xi in x]

err1 = [abs(y1[i] - y_an[i]) for i in range(n)]
err2 = [abs(y2[i] - y_an[i]) for i in range(n)]

print(f"{'x':>6} {'1-й порядок':>12} {'2-й порядок':>12} {'Аналит.':>12} {'Ошибка 1':>12} {'Ошибка 2':>12}")
for i in range(n):
    print(f"{x[i]:6.1f} {y1[i]:12.6f} {y2[i]:12.6f} {y_an[i]:12.6f} {err1[i]:12.2e} {err2[i]:12.2e}")

print(f"\nМакс. погрешность (1-й порядок): {max(err1):.2e}")
print(f"Макс. погрешность (2-й порядок): {max(err2):.2e}")

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x, y_an, 'k-', label='Аналитическое', linewidth=2)
plt.plot(x, y1, 'ro--', label='1-й порядок', markersize=5)
plt.plot(x, y2, 'bs--', label='2-й порядок', markersize=5)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.title('Решения')

plt.subplot(1, 2, 2)
plt.semilogy(x, err1, 'ro-', label='1-й порядок')
plt.semilogy(x, err2, 'bs-', label='2-й порядок')
plt.xlabel('x')
plt.ylabel('Абсолютная погрешность')
plt.legend()
plt.grid(True)
plt.title('Погрешности')

plt.tight_layout()
plt.show()
