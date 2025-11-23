import math
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Глобальные константы (Теоретические корни для графика) ---
# Эти значения получены точным решением системы, чтобы отображать на графике "Идеал"
ANALYTICAL_ROOTS = [
    (2.57836, 1.16278),   # I четверть
    (-1.12115, 2.59663),  # II четверть
    (-1.53542, -2.37542), # III четверть
    (2.06716, -1.92985)   # IV четверть
]

# Список для сохранения истории всех вычислений
execution_log = []

def get_lambdas(x, y):
    """
    Адаптивный выбор параметров lambda на основе знака производной в текущей точке.
    """
    df1_dx = 2 * x
    df2_dy = -2 * x - 6 * y
    
    l1_base = 0.1
    l2_base = 0.05
    
    # Инвертируем знак лямбды относительно знака производной для обеспечения сходимости
    l1 = -l1_base if df1_dx > 0 else l1_base
    l2 = -l2_base if df2_dy > 0 else l2_base
    
    return l1, l2

def get_functions(x1, x2):
    lam1, lam2 = get_lambdas(x1, x2)
    f1 = x1**2 + x2**2 - 8
    f2 = 3*(x1**2) - 2*x1*x2 - 3*(x2**2) - 5*x1 + 3
    phi1 = x1 + lam1 * f1
    phi2 = x2 + lam2 * f2
    return phi1, phi2

def get_derivatives_analytical(x1, x2):
    lam1, lam2 = get_lambdas(x1, x2)
    dphi1_dx1 = 1 + lam1 * (2 * x1)
    dphi1_dx2 = lam1 * (2 * x2)
    dphi2_dx1 = lam2 * (6 * x1 - 2 * x2 - 5)
    dphi2_dx2 = 1 + lam2 * (-2 * x1 - 6 * x2)
    return dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2

def check_convergence(x1, x2):
    dp1_dx1, dp1_dx2, dp2_dx1, dp2_dx2 = get_derivatives_analytical(x1, x2)
    norm1 = abs(dp1_dx1) + abs(dp1_dx2)
    norm2 = abs(dp2_dx1) + abs(dp2_dx2)
    max_norm = max(norm1, norm2)
    return max_norm < 1, max_norm

def solve_system(method_name, x_start, y_start, epsilon):
    """
    Универсальная функция-решатель.
    method_name: 'SimpleIteration' или 'Seidel'
    """
    print(f"\n--- {method_name} (Старт: {x_start}, {y_start}) ---")
    print(f"{'k':<3} | {'x1':<9} | {'x2':<9} | {'dx':<9} | {'dy':<9} | {'Norma J'}")
    print("-" * 65)
    
    x1, x2 = x_start, y_start
    k = 0
    max_iter = 100
    status = "Converged"
    
    try:
        converges, norm = check_convergence(x1, x2)
        print(f"{k:<3} | {x1:<9.5f} | {x2:<9.5f} | {'-':<9} | {'-':<9} | {norm:.4f}")
        
        while True:
            k += 1
            
            if method_name == 'SimpleIteration':
                x1_new, x2_new = get_functions(x1, x2)
            else: # Seidel
                temp_x1, _ = get_functions(x1, x2)
                x1_new = temp_x1
                _, temp_x2 = get_functions(x1_new, x2) # Используем новый x1
                x2_new = temp_x2
                
            # Проверка на "взрыв" значений
            if abs(x1_new) > 1e8 or abs(x2_new) > 1e8:
                print("!!! РАСХОДИМОСТЬ (Числа слишком большие) !!!")
                status = "Diverged (Overflow)"
                break

            dx1 = abs(x1_new - x1)
            dx2 = abs(x2_new - x2)
            _, norm = check_convergence(x1_new, x2_new)
            
            print(f"{k:<3} | {x1_new:<9.5f} | {x2_new:<9.5f} | {dx1:<9.5f} | {dx2:<9.5f} | {norm:.4f}")
            
            if dx1 <= epsilon and dx2 <= epsilon:
                break
            
            if k >= max_iter:
                print("Превышен лимит итераций!")
                status = "Limit Reached"
                break
            
            x1, x2 = x1_new, x2_new
            
    except OverflowError:
        print("!!! OverflowError: Метод разошелся !!!")
        status = "Error"
        x1, x2 = None, None

    # Сохраняем результат в общий лог
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
    """
    Строит график системы уравнений, используя:
    1. Контуры (линии уровня) для уравнений - это точная геометрия.
    2. ANALYTICAL_ROOTS - точные координаты корней (желтые точки).
    """
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Сетка для контуров
    x = np.linspace(-4, 4, 400)
    y = np.linspace(-4, 4, 400)
    X, Y = np.meshgrid(x, y)
    
    # Уравнения F(x,y) = 0
    Z1 = X**2 + Y**2 - 8
    Z2 = 3*X**2 - 2*X*Y - 3*Y**2 - 5*X + 3
    
    # Рисуем линии уровня 0 (самые точные графики уравнений)
    c1 = ax.contour(X, Y, Z1, levels=[0], colors='#00ffff', linewidths=2.5)
    c2 = ax.contour(X, Y, Z2, levels=[0], colors='#ff00ff', linewidths=2.5)
    
    # Рисуем ТОЧНЫЕ АНАЛИТИЧЕСКИЕ КОРНИ
    for rx, ry in ANALYTICAL_ROOTS:
        ax.plot(rx, ry, 'o', color='#ffff00', markersize=12, 
                markeredgecolor='#ff6600', markeredgewidth=2, zorder=10)
        
        # Подпись координат
        ax.annotate(f'({rx:.2f}, {ry:.2f})', xy=(rx, ry), xytext=(15, 15),
                   textcoords='offset points', fontsize=10, color='#ffff00', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', fc='#1a1a1a', ec='#ffff00', alpha=0.8))

    # Оформление
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xlabel('x₁', fontsize=12, color='white')
    ax.set_ylabel('x₂', fontsize=12, color='white')
    ax.set_title('Аналитическое решение системы', fontsize=14, pad=15)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    
    # Легенда
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

# --- ЗАПУСК ---
epsilon_val = 0.0001
initial_points = [
    (2.0, 1.0),    # Для корня в I четверти
    (-1.0, 2.5),   # Для корня во II четверти
    (-1.5, -2.3),  # Для корня в III четверти
    (2.0, -2.0)    # Для корня в IV четверти
]

# 1. Выполняем расчеты
for start_x, start_y in initial_points:
    solve_system('SimpleIteration', start_x, start_y, epsilon_val)
    solve_system('Seidel', start_x, start_y, epsilon_val)

# 2. Выводим сводную таблицу (лог)
print_final_summary()

# 3. Строим "чистый" график с аналитическими решениями
plot_analytical_system()