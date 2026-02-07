import matplotlib.pyplot as plt

def analytical_solution(x):
    return (x * (x - 2)) ** 0.5 + x ** 0.5

def p(x):
    return -x * (x - 4) / (2 * x**2 * (x - 2))

def q(x):
    return (x - 3) / (2 * x**2 * (x - 2))

def gauss_solve(A, b):
    n = len(b)
    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > abs(A[max_row][i]):
                max_row = k
        A[i], A[max_row] = A[max_row], A[i]
        b[i], b[max_row] = b[max_row], b[i]
        
        for k in range(i + 1, n):
            factor = A[k][i] / A[i][i]
            b[k] -= factor * b[i]
            for j in range(i, n):
                A[k][j] -= factor * A[i][j]
    
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= A[i][j] * x[j]
        x[i] /= A[i][i]
    return x

def solve_boundary_problem(x_start, x_end, h, order=1):
    x_vals = []
    x = x_start
    while x <= x_end + 1e-9:
        x_vals.append(x)
        x += h
    
    n = len(x_vals)
    A = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    
    for i in range(1, n - 1):
        xi = x_vals[i]
        pi = p(xi)
        qi = q(xi)
        
        A[i][i - 1] = 1 / h**2 - pi / (2 * h)
        A[i][i] = -2 / h**2 + qi
        A[i][i + 1] = 1 / h**2 + pi / (2 * h)
        b[i] = 0
    
    if order == 1:
        A[0][0] = 3 + 1 / h
        A[0][1] = -1 / h
        b[0] = 9.755
    else:
        A[0][0] = 3 - 3 / (2 * h)
        A[0][1] = 4 / (2 * h)
        A[0][2] = -1 / (2 * h)
        b[0] = 9.755
    
    A[-1][-1] = 1
    b[-1] = 10.937
    
    y = gauss_solve(A, b)
    return x_vals, y

x_start, x_end, h = 2.5, 9, 0.5

x1, y1 = solve_boundary_problem(x_start, x_end, h, order=1)
x2, y2 = solve_boundary_problem(x_start, x_end, h, order=2)

y_analytical = [analytical_solution(x) for x in x1]

error1 = [abs(y1[i] - y_analytical[i]) for i in range(len(y1))]
error2 = [abs(y2[i] - y_analytical[i]) for i in range(len(y2))]

max_error1 = max(error1)
max_error2 = max(error2)

print("Порядок 1: макс. погрешность =", max_error1)
print("Порядок 2: макс. погрешность =", max_error2)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x1, y_analytical, 'k-', label='Аналитическое', linewidth=2)
plt.plot(x1, y1, 'ro--', label='Порядок 1', markersize=5)
plt.plot(x2, y2, 'bs--', label='Порядок 2', markersize=5)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.title('Решения')

plt.subplot(1, 2, 2)
plt.plot(x1, error1, 'ro-', label='Порядок 1', markersize=5)
plt.plot(x2, error2, 'bs-', label='Порядок 2', markersize=5)
plt.xlabel('x')
plt.ylabel('Погрешность')
plt.legend()
plt.grid(True)
plt.title('Погрешности')

plt.tight_layout()
plt.show()
