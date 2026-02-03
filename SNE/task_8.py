import math
import matplotlib.pyplot as plt
import numpy as np

eps = 0.00001

ANALYTICAL_ROOTS = [
    (2.5784, 1.1627),   # I четверть
    (-1.2091, 2.5569),  # II четверть
    (-1.8163, -2.1681), # III четверть
    (1.9470, -2.0517)   # IV четверть
]

# Список для сохранения истории всех вычислений
execution_log = []

def f1(x1, x2):
    return x1**2 + x2**2 - 8

def f2(x1, x2):
    return 3*x1**2 - 2*x1*x2 - 3*x2**2 - 5*x1 + 3

def df1_dx1(x1, x2):
    return 2 * x1

def df1_dx2(x1, x2):
    return 2 * x2

def df2_dx1(x1, x2):
    return 6*x1 - 2*x2 - 5

def df2_dx2(x1, x2):
    return -2*x1 - 6*x2

def get_lambdas(x, y):
    l1_base = 0.1
    l2_base = 0.05
    
    l1 = -l1_base if df1_dx1(x, y) > 0 else l1_base
    l2 = -l2_base if df2_dx2(x, y) > 0 else l2_base
    
    return l1, l2

def phi1(x1, x2):
    lam1, _ = get_lambdas(x1, x2)
    return x1 + lam1 * f1(x1, x2)

def phi2(x1, x2):
    _, lam2 = get_lambdas(x1, x2)
    return x2 + lam2 * f2(x1, x2)

def dphi1_dx1(x1, x2):
    lam1, _ = get_lambdas(x1, x2)
    return 1 + lam1 * df1_dx1(x1, x2)

def dphi1_dx2(x1, x2):
    lam1, _ = get_lambdas(x1, x2)
    return lam1 * df1_dx2(x1, x2)

def dphi2_dx1(x1, x2):
    _, lam2 = get_lambdas(x1, x2)
    return lam2 * df2_dx1(x1, x2)

def dphi2_dx2(x1, x2):
    _, lam2 = get_lambdas(x1, x2)
    return 1 + lam2 * df2_dx2(x1, x2)

def jacobian(x1, x2):
    return df1_dx1(x1, x2), df1_dx2(x1, x2), df2_dx1(x1, x2), df2_dx2(x1, x2)

def check_convergence(x1, x2):
    norm_row1 = abs(dphi1_dx1(x1, x2)) + abs(dphi1_dx2(x1, x2))
    norm_row2 = abs(dphi2_dx1(x1, x2)) + abs(dphi2_dx2(x1, x2))
    max_norm_rows = max(norm_row1, norm_row2)
    
    norm_col1 = abs(dphi1_dx1(x1, x2)) + abs(dphi2_dx1(x1, x2))
    norm_col2 = abs(dphi1_dx2(x1, x2)) + abs(dphi2_dx2(x1, x2))
    max_norm_cols = max(norm_col1, norm_col2)
    
    max_norm = max(max_norm_rows, max_norm_cols)
    return max_norm < 1, max_norm

def simple_iteration(x1, x2, epsilon, max_iter=100):
    k = 0
    iterations = [(0, x1, x2, None, None)]
    while k < max_iter:
        k += 1
        x1_new, x2_new = phi1(x1, x2), phi2(x1, x2)
        
        if abs(x1_new) > 1e8 or abs(x2_new) > 1e8:
            return None, None, k, "Diverged (Overflow)", iterations
        
        dx1, dx2 = abs(x1_new - x1), abs(x2_new - x2)
        iterations.append((k, x1_new, x2_new, dx1, dx2))
        
        if dx1 <= eps and dx2 <= eps:
            return x1_new, x2_new, k, "Converged", iterations
        
        x1, x2 = x1_new, x2_new
    
    return x1, x2, k, "Limit Reached", iterations

def seidel(x1, x2, epsilon, max_iter=100):
    k = 0
    iterations = [(0, x1, x2, None, None)]
    while k < max_iter:
        k += 1
        x1_new = phi1(x1, x2)
        x2_new = phi2(x1_new, x2)
        
        if abs(x1_new) > 1e8 or abs(x2_new) > 1e8:
            return None, None, k, "Diverged (Overflow)", iterations
        
        dx1, dx2 = abs(x1_new - x1), abs(x2_new - x2)
        iterations.append((k, x1_new, x2_new, dx1, dx2))
        
        if dx1 <= eps and dx2 <= eps:
            return x1_new, x2_new, k, "Converged", iterations
        
        x1, x2 = x1_new, x2_new
    
    return x1, x2, k, "Limit Reached", iterations

def newton(x1, x2, epsilon, max_iter=100):
    k = 0
    iterations = [(0, x1, x2, None, None)]
    while k < max_iter:
        k += 1
        f1_val = f1(x1, x2)
        f2_val = f2(x1, x2)
        
        j11, j12, j21, j22 = jacobian(x1, x2)
        
        det = j11*j22 - j12*j21
        if abs(det) < 1e-10:
            return None, None, k, "Singular Jacobian", iterations
        
        dx1 = (j22*f1_val - j12*f2_val) / det
        dx2 = (j11*f2_val - j21*f1_val) / det
        
        x1_new, x2_new = x1 - dx1, x2 - dx2
        
        if abs(x1_new) > 1e8 or abs(x2_new) > 1e8:
            return None, None, k, "Diverged (Overflow)", iterations
        
        iterations.append((k, x1_new, x2_new, abs(dx1), abs(dx2)))
        
        if abs(dx1) <= eps and abs(dx2) <= eps:
            return x1_new, x2_new, k, "Converged", iterations
        
        x1, x2 = x1_new, x2_new
    
    return x1, x2, k, "Limit Reached", iterations


def solve_system(method_name, x_start, y_start, epsilon):
    print(f"\n--- {method_name} ---")
    print(f"{'k':<3} | {'x1':<9} | {'x2':<9} | {'dx':<9} | {'dy':<9}")
    print("-" * 50)
    
    x1, x2 = x_start, y_start
    
    try:
        if method_name == 'SimpleIteration':
            x1, x2, k, status, iterations = simple_iteration(x1, x2, epsilon)
        elif method_name == 'Seidel':
            x1, x2, k, status, iterations = seidel(x1, x2, epsilon)
        else:
            x1, x2, k, status, iterations = newton(x1, x2, epsilon)
        
        for iter_k, iter_x1, iter_x2, iter_dx1, iter_dx2 in iterations:
            if iter_dx1 is None:
                print(f"{iter_k:<3} | {iter_x1:<9.5f} | {iter_x2:<9.5f} | {'-':<9} | {'-':<9}")
            else:
                print(f"{iter_k:<3} | {iter_x1:<9.5f} | {iter_x2:<9.5f} | {iter_dx1:<9.5f} | {iter_dx2:<9.5f}")
        
        if x1 is None:
            print(f"Статус: {status}")
            
    except OverflowError:
        status = "Error"
        x1, x2 = None, None
        k = 0

    execution_log.append({
        'method': method_name,
        'start': (x_start, y_start),
        'end': (x1, x2) if x1 is not None else (0,0),
        'iters': k,
        'status': status
    })
    
    return x1, x2

def print_final_summary():
    print("\n" + "="*80)
    print(f"{'ИТОГОВЫЙ ОТЧЕТ':^80}")
    print("="*80)
    print(f"{'Метод':<20} | {'Старт':<15} | {'Корень (x, y)':<20} | {'Итер.':<5} | {'Статус'}")
    print("-" * 80)
    
    for entry in execution_log:
        start_s = f"({entry['start'][0]:.1f}, {entry['start'][1]:.1f})"
        if entry['end'] and entry['end'][0] is not None:
            end_s = f"({entry['end'][0]:.4f}, {entry['end'][1]:.4f})"
        else:
            end_s = "None"
            
        print(f"{entry['method']:<20} | {start_s:<15} | {end_s:<20} | {entry['iters']:<5} | {entry['status']}")
    print("="*80)

def plot_analytical_system():
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    x = np.linspace(-4, 4, 400)
    y = np.linspace(-4, 4, 400)
    X, Y = np.meshgrid(x, y)
    
    Z1 = X**2 + Y**2 - 8
    Z2 = 3*X**2 - 2*X*Y - 3*Y**2 - 5*X + 3
    
    c1 = ax.contour(X, Y, Z1, levels=[0], colors='#00ffff', linewidths=2.5)
    c2 = ax.contour(X, Y, Z2, levels=[0], colors='#ff00ff', linewidths=2.5)
    
    for rx, ry in ANALYTICAL_ROOTS:
        ax.plot(rx, ry, 'o', color='#ffff00', markersize=12, 
                markeredgecolor='#ff6600', markeredgewidth=2, zorder=10)
        
        ax.annotate(f'({rx:.2f}, {ry:.2f})', xy=(rx, ry), xytext=(15, 15),
                   textcoords='offset points', fontsize=10, color='#ffff00', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', fc='#1a1a1a', ec='#ffff00', alpha=0.8))

    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xlabel('x₁', fontsize=12, color='white')
    ax.set_ylabel('x₂', fontsize=12, color='white')
    ax.set_title('Аналитическое решение системы', fontsize=14, pad=15)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='#00ffff', lw=2, label='x² + y² - 8 = 0'),
        Line2D([0], [0], color='#ff00ff', lw=2, label='3x² - 2xy - 3y² - 5x + 3 = 0'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ffff00', 
               markeredgecolor='#ff6600', markersize=10, label='Точные корни')
    ]
    ax.legend(handles=legend_elements, loc='upper right', facecolor='#1a1a1a')
    
    plt.tight_layout()
    plt.show()

initial_points = [
    (2.0, 1.0),
    (-1.0, 2.5),
    (-1.5, -2.3),
    (2.0, -2.0)
]

# 1. Выполняем расчеты
for i, (start_x, start_y) in enumerate(initial_points, 1):
    print(f"\n{'='*60}")
    print(f"КОРЕНЬ #{i}: Начальное приближение ({start_x}, {start_y})")
    converges, norm = check_convergence(start_x, start_y)
    print(f"Проверка сходимости: ||J|| = {norm:.4f} {'< 1 ✓' if converges else '>= 1 ✗'}")
    print(f"{'='*60}")
    
    solve_system('SimpleIteration', start_x, start_y, eps)
    solve_system('Seidel', start_x, start_y, eps)
    solve_system('Newton', start_x, start_y, eps)

# 2. Выводим сводную таблицу (лог)
print_final_summary()

# 3. Строим "чистый" график с аналитическими решениями
plot_analytical_system()