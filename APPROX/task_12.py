import matplotlib.pyplot as plt

X = [4.71, 3.63, 2.55, 1.47, 0.39, 0.69, 1.77, 2.85, 3.93, 5.01, 6.09]
Y = [2.9853, 1.7984, 1.0358, 0.6987, 0.7116, 1.3214, 1.5742, 1.6371, 1.3936, 1.0309, 0.9754]
x_star = -0.926

def gauss(A, b):
    n = len(b)
    for i in range(n):
        for j in range(i+1, n):
            m = A[j][i] / A[i][i]
            for k in range(n):
                A[j][k] -= m * A[i][k]
            b[j] -= m * b[i]
    x = [0] * n
    for i in range(n-1, -1, -1):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i+1, n))) / A[i][i]
    return x

def least_squares(X, Y, degree):
    n = len(X)
    m = degree + 1
    A = [[sum(X[k]**(i+j) for k in range(n)) for j in range(m)] for i in range(m)]
    b = [sum(Y[k] * X[k]**i for k in range(n)) for i in range(m)]
    return gauss(A, b)

def poly(x, c):
    return sum(c[i] * x**i for i in range(len(c)))

def sse(X, Y, c):
    return sum((Y[i] - poly(X[i], c))**2 for i in range(len(X)))

c1 = least_squares(X, Y, 1)
c2 = least_squares(X, Y, 2)
c3 = least_squares(X, Y, 3)

e1, e2, e3 = sse(X, Y, c1), sse(X, Y, c2), sse(X, Y, c3)

print(f"P1(x) = {c1[0]:.4f} + {c1[1]:.4f}*x")
print(f"Ошибка: {e1:.6f}, P1({x_star}) = {poly(x_star, c1):.6f}\n")

print(f"P2(x) = {c2[0]:.4f} + {c2[1]:.4f}*x + {c2[2]:.4f}*x^2")
print(f"Ошибка: {e2:.6f}, P2({x_star}) = {poly(x_star, c2):.6f}\n")

print(f"P3(x) = {c3[0]:.4f} + {c3[1]:.4f}*x + {c3[2]:.4f}*x^2 + {c3[3]:.4f}*x^3")
print(f"Ошибка: {e3:.6f}, P3({x_star}) = {poly(x_star, c3):.6f}\n")

x_plot = [min(min(X), x_star) - 0.5 + i * 0.01 for i in range(750)]

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for ax, c, deg, col in [(axes[0,0], c1, 1, 'b'), (axes[0,1], c2, 2, 'g'), (axes[1,0], c3, 3, 'r')]:
    y_plot = [poly(x, c) for x in x_plot]
    ax.plot(X, Y, 'ko', markersize=6)
    ax.plot(x_plot, y_plot, col+'-')
    ax.plot(x_star, poly(x_star, c), col+'*', markersize=12)
    ax.axvline(x=x_star, color='gray', linestyle='--', alpha=0.5)
    ax.set_title(f'Многочлен степени {deg}')
    ax.grid(True, alpha=0.3)

axes[1,1].bar(['Степень 1', 'Степень 2', 'Степень 3'], [e1, e2, e3], color=['b', 'g', 'r'])
axes[1,1].set_title('Сравнение точности')
axes[1,1].set_ylabel('Сумма квадратов ошибок')
axes[1,1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

