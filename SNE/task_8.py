import math
import matplotlib.pyplot as plt
import numpy as np

# Система уравнений:
# F1: x1² + x2² - 8 = 0
# F2: 3x1² - 2x1x2 - 3x2² - 5x1 + 3 = 0

def F1(x, y):
    return x*x + y*y - 8

def F2(x, y):
    return 3*x*x - 2*x*y - 3*y*y - 5*x + 3

def dF1_dx(x, y):
    return 2*x

def dF1_dy(x, y):
    return 2*y

def dF2_dx(x, y):
    return 6*x - 2*y - 5

def dF2_dy(x, y):
    return -2*x - 6*y

def jacobian(x, y):
    return [[dF1_dx(x, y), dF1_dy(x, y)],
            [dF2_dx(x, y), dF2_dy(x, y)]]

def det2x2(m):
    return m[0][0]*m[1][1] - m[0][1]*m[1][0]

def inverse2x2(m):
    d = det2x2(m)
    if abs(d) < 1e-10:
        return None
    return [[m[1][1]/d, -m[0][1]/d],
            [-m[1][0]/d, m[0][0]/d]]

def mul_mat_vec(m, v):
    return [m[0][0]*v[0] + m[0][1]*v[1],
            m[1][0]*v[0] + m[1][1]*v[1]]

def mul_mat_mat(a, b):
    return [[a[0][0]*b[0][0] + a[0][1]*b[1][0], a[0][0]*b[0][1] + a[0][1]*b[1][1]],
            [a[1][0]*b[0][0] + a[1][1]*b[1][0], a[1][0]*b[0][1] + a[1][1]*b[1][1]]]

def vec_norm(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1])

def newton_system(x0, y0, eps, max_iter):
    x, y = x0, y0
    
    for i in range(max_iter):
        J = jacobian(x, y)
        Jinv = inverse2x2(J)
        if Jinv is None:
            return None, None, 0, False, "Матрица Якоби вырожденная"
        
        f1, f2 = F1(x, y), F2(x, y)
        delta = mul_mat_vec(Jinv, [-f1, -f2])
        x += delta[0]
        y += delta[1]
        
        if vec_norm(delta) < eps:
            return x, y, i+1, True, "Метод сходится"
    
    return x, y, max_iter, False, "Достигнуто максимальное количество итераций"

def simple_iteration_system(x0, y0, eps, max_iter):
    x, y = x0, y0
    tau = 0.5
    
    # Проверка условия сходимости
    sup_norm = 0.0
    samples = 5
    radius = 0.5
    
    for i in range(samples):
        for j in range(samples):
            xi = x0 + (-radius + (2*radius)*i/(samples-1))
            yi = y0 + (-radius + (2*radius)*j/(samples-1))
            
            J = jacobian(xi, yi)
            Jinv = inverse2x2(J)
            if Jinv:
                # Phi'(x) = I - tau*Jinv*J = I - tau*I = (1-tau)*I
                # Но на самом деле Phi(x) = x - tau*Jinv*F, поэтому Phi'(x) = I - tau*Jinv*J
                JinvJ = mul_mat_mat(Jinv, J)
                phi_deriv = [[1 - tau*JinvJ[0][0], -tau*JinvJ[0][1]],
                            [-tau*JinvJ[1][0], 1 - tau*JinvJ[1][1]]]
                row1_norm = abs(phi_deriv[0][0]) + abs(phi_deriv[0][1])
                row2_norm = abs(phi_deriv[1][0]) + abs(phi_deriv[1][1])
                max_row_norm = max(row1_norm, row2_norm)
                if max_row_norm > sup_norm:
                    sup_norm = max_row_norm
    
    for k in range(max_iter):
        J = jacobian(x, y)
        Jinv = inverse2x2(J)
        if Jinv is None:
            return None, None, 0, False, "Матрица Якоби вырожденная"
        
        f = [F1(x, y), F2(x, y)]
        delta = mul_mat_vec(Jinv, [-f[0], -f[1]])
        
        x_new = x + tau * delta[0]
        y_new = y + tau * delta[1]
        
        if max(abs(x_new - x), abs(y_new - y)) < eps:
            return x_new, y_new, k+1, True, f"sup||Φ'|| ~= {sup_norm}, решение найдено за {k+1} итераций"
        
        x, y = x_new, y_new
    
    return x, y, max_iter, False, f"sup||Φ'|| ~= {sup_norm}, не сошёлся за {max_iter} итераций"

def seidel_system(x0, y0, eps, max_iter):
    x, y = x0, y0
    tau = 0.5
    
    # Проверка условия сходимости
    sup_norm = 0.0
    samples = 5
    radius = 0.5
    
    for i in range(samples):
        for j in range(samples):
            xi = x0 + (-radius + (2*radius)*i/(samples-1))
            yi = y0 + (-radius + (2*radius)*j/(samples-1))
            
            J = jacobian(xi, yi)
            Jinv = inverse2x2(J)
            if Jinv:
                JinvJ = mul_mat_mat(Jinv, J)
                phi_deriv = [[1 - tau*JinvJ[0][0], -tau*JinvJ[0][1]],
                            [-tau*JinvJ[1][0], 1 - tau*JinvJ[1][1]]]
                row1_norm = abs(phi_deriv[0][0]) + abs(phi_deriv[0][1])
                row2_norm = abs(phi_deriv[1][0]) + abs(phi_deriv[1][1])
                max_row_norm = max(row1_norm, row2_norm)
                if max_row_norm > sup_norm:
                    sup_norm = max_row_norm
    
    for k in range(max_iter):
        J = jacobian(x, y)
        Jinv = inverse2x2(J)
        if Jinv is None:
            return None, None, 0, False, "Матрица Якоби вырожденная"
        
        f1 = F1(x, y)
        delta1 = mul_mat_vec(Jinv, [-f1, 0])
        x_new = x + tau * delta1[0]
        
        f2 = F2(x_new, y)
        J_new = jacobian(x_new, y)
        Jinv_new = inverse2x2(J_new)
        if Jinv_new is None:
            return None, None, 0, False, "Матрица Якоби вырожденная"
        
        delta2 = mul_mat_vec(Jinv_new, [0, -f2])
        y_new = y + tau * delta2[1]
        
        if max(abs(x_new - x), abs(y_new - y)) < eps:
            return x_new, y_new, k+1, True, f"sup||Φ'|| ~= {sup_norm:.3f}, решение найдено за {k+1} итераций"
        
        x, y = x_new, y_new
    
    return x, y, max_iter, False, f"sup||Φ'|| ~= {sup_norm:.3f}, не сошёлся за {max_iter} итераций"

def plot_system(solutions):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 10))
    
    x = np.linspace(-4, 4, 1000)
    y = np.linspace(-4, 4, 1000)
    X, Y = np.meshgrid(x, y)
    
    Z1 = X**2 + Y**2 - 8
    Z2 = 3*X**2 - 2*X*Y - 3*Y**2 - 5*X + 3
    
    c1 = ax.contour(X, Y, Z1, levels=[0], colors='#00ffff', linewidths=3, linestyles='-')
    c2 = ax.contour(X, Y, Z2, levels=[0], colors='#ff00ff', linewidths=3, linestyles='-')
    
    for sol_x, sol_y in solutions:
        ax.plot(sol_x, sol_y, 'o', color='#ffff00', markersize=15, markeredgecolor='#ff6600', markeredgewidth=3, zorder=5)
        ax.annotate(f'({sol_x:.2f}, {sol_y:.2f})', xy=(sol_x, sol_y), xytext=(10, 10),
                   textcoords='offset points', fontsize=11, color='#ffff00',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='#1a1a1a', edgecolor='#ffff00', linewidth=2))
    
    ax.grid(True, alpha=0.3, linestyle='--', color='#404040')
    ax.axhline(y=0, color='#808080', linewidth=1.5, alpha=0.5)
    ax.axvline(x=0, color='#808080', linewidth=1.5, alpha=0.5)
    ax.set_xlabel('x₁', fontsize=14, color='#00ffff', fontweight='bold')
    ax.set_ylabel('x₂', fontsize=14, color='#ff00ff', fontweight='bold')
    ax.set_title('Система нелинейных уравнений', fontsize=16, color='#ffffff', fontweight='bold', pad=20)
    
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='#00ffff', linewidth=3, label='F₁ = 0'),
                      Line2D([0], [0], color='#ff00ff', linewidth=3, label='F₂ = 0'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='#ffff00', 
                            markeredgecolor='#ff6600', markeredgewidth=3, markersize=10, label='Решения', linestyle='None')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.9, facecolor='#1a1a1a', edgecolor='#ffffff')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    eps = 0.0001
    max_iter = 1000

    initial_values = [(2.0, 1.0), (-1.0, 2.5), (-1.5, -2.0), (2.0, -2.0)]
    solutions = []
    
    print("Система уравнений:")
    print("F₁: x₁² + x₂² - 8 = 0")
    print("F₂: 3x₁² - 2x₁x₂ - 3x₂² - 5x₁ + 3 = 0\n")

    for value in initial_values:
        x0, y0 = value
        print(f"=== Начальное приближение x₀={x0}, y₀={y0} ===")
        print("=== Метод Ньютона ===")
        x, y, iters, conv, msg = newton_system(x0, y0, eps, max_iter)
        print(f"x = {x:.6f}, y = {y:.6f}")
        print(f"Итераций: {iters}")
        print(f"{msg}\n")
        
        if conv and (x, y) not in [(s[0], s[1]) for s in solutions]:
            solutions.append((x, y))

        print("=== Метод простой итерации ===")
        x, y, iters, conv, msg = simple_iteration_system(x0, y0, eps, max_iter)
        print(f"x = {x:.6f}, y = {y:.6f}")
        print(f"Итераций: {iters}")
        print(f"{msg}\n")

        print("=== Метод Зейделя ===")
        x, y, iters, conv, msg = seidel_system(x0, y0, eps, max_iter)
        print(f"x = {x:.6f}, y = {y:.6f}")
        print(f"Итераций: {iters}")
        print(f"{msg}")

        print("============\n\n\n")
    
    plot_system(solutions)
