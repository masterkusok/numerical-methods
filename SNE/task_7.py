import math
import matplotlib.pyplot as plt

eps = 0.0001

# Уравнение: e^(x^2) + 3^(-x) - 5x^2 - 1 = 0
def f(x):
    """Функция уравнения"""
    try:
        return math.exp(x**2) + 3**(-x) - 5*x**2 - 1
    except:
        return float('inf')

def df(x):
    """Первая производная"""
    try:
        return 2*x*math.exp(x**2) - 3**(-x)*math.log(3) - 10*x
    except:
        return float('inf')

def d2f(x):
    """Вторая производная"""
    try:
        return (2 + 4*x**2)*math.exp(x**2) + 3**(-x)*(math.log(3))**2 - 10
    except:
        return float('inf')

def phi_lambda(x, lam):
    sign = 1 if df(x) > 0 else -1
    return x - sign * lam * f(x)

def phi_derivative_lambda(x, lam):
    sign = 1 if df(x) > 0 else -1
    return 1 - sign * lam * df(x)

def simple_iteration(x0, a, b, lam, max_iter=1000):
    phi_a = abs(phi_derivative_lambda(a, lam))
    phi_b = abs(phi_derivative_lambda(b, lam))
    phi_x0 = abs(phi_derivative_lambda(x0, lam))

    if phi_a < 1.0 and phi_b < 1.0 and phi_x0 < 1.0:
        print("Выполнено условие сходимости!")
    else:
        print("Условие сходимости НЕ выполнено!")
        print(f"   |φ'({a})| = {phi_a:.6f}, |φ'({b})| = {phi_b:.6f}")
        return 0, 0, []

    iterations = [x0]
    x = x0

    for i in range(max_iter):
        x_new  = phi_lambda(x, lam)
        iterations.append(x_new)

        if abs(x_new - x) < eps:
            return x_new, i + 1, iterations
        x = x_new

    return x, max_iter, iterations

# 2. Метод Ньютона
def newton_method(a, b, max_iter=1000):
    x0 = a if f(a) * d2f(a) > 0 else b
    
    iterations = []
    x = x0
    iterations.append(x)
    
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)
        
        if abs(dfx) < eps:
            print("   Производная близка к нулю")
            break
            
        x_new = x - fx / dfx
        iterations.append(x_new)
        
        if abs(x_new - x) < eps:
            return x_new, i + 1, iterations
        x = x_new
    
    return x, max_iter, iterations

# 3. Метод секущих
def secant_method(a, b, max_iter=1000):
    x0 = a if f(a) * d2f(a) > 0 else b
    
    iterations = []
    iterations.append(x0)
    
    fx0 = f(x0)
    dfx0 = df(x0)
    
    if dfx0 == 0:
        print("Ошибка: Производная в начальной точке равна нулю.")
        return None, 0, iterations

    x1 = x0 - fx0 / dfx0
    iterations.append(x1)
    
    for i in range(max_iter):
        fx1 = f(x1)
        fx0 = f(x0)
        
        if abs(fx1 - fx0) < 1e-12:
            print("   Деление на ноль в методе секущих (значения функции совпали)")
            break
        
        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        iterations.append(x_new)
        
        if abs(x_new - x1) < eps:
            return x_new, i + 1, iterations
        
        x0 = x1
        x1 = x_new
    
    return x1, max_iter, iterations

# 4. Метод хорд
def chord_method(a, b, max_iter=1000):
    """Метод хорд (метод ложного положения)"""
    iterations = []
    
    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        print("   ОШИБКА: f(a)*f(b) > 0 - корня на отрезке нет!")
        return None, 0, []
    
    # Выбор начального приближения и неподвижной точки
    if f(a) * d2f(a) > 0:
        x = b  # Начальное приближение
        z = a  # Неподвижная точка
    else:
        x = a  # Начальное приближение
        z = b  # Неподвижная точка

    fz = f(z)
    iterations.append(x)
    
    for i in range(max_iter):
        fx = f(x)
        
        if abs(fx - fz) < 1e-12:
            break
            
        x_new = x - fx * (x - z) / (fx - fz)
        iterations.append(x_new)
        
        if abs(x_new - x) < eps:
            return x_new, i + 1, iterations
        x = x_new
    
    return x, max_iter, iterations

# 5. Метод дихотомии (деления пополам)
def bisection_method(a, b, max_iter=1000):
    """Метод дихотомии"""
    iterations = []
    
    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        print("   ОШИБКА: f(a)*f(b) > 0 - корня на отрезке нет!")
        return None, 0, []
    
    for i in range(max_iter):
        c = (a + b) / 2
        iterations.append(c)
        fc = f(c)
        
        if abs(fc) < eps or abs(b - a) < 2*eps:
            return c, i + 1, iterations
        
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    
    return (a + b) / 2, max_iter, iterations

def plot_function():
    """Построение графика функции"""
    x_vals = [i * 0.01 for i in range(-200, 200)]
    y_vals = []
    
    for x in x_vals:
        try:
            y = f(x)
            if abs(y) < 50:
                y_vals.append(y)
            else:
                y_vals.append(None)
        except:
            y_vals.append(None)
    
    plt.figure(figsize=(14, 5))
    plt.plot(x_vals, y_vals, 'b-', linewidth=2, label='$e^{x^2} + 3^{-x} - 5x^2 - 1$')
    plt.axhline(y=0, color='r', linestyle='--', linewidth=1, label='y = 0')
    plt.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    plt.grid(True, alpha=0.3)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('f(x)', fontsize=12)
    plt.title('График функции $f(x) = e^{x^2} + 3^{-x} - 5x^2 - 1$', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.ylim(-15, 15)
    
    plt.plot([-0.85, 0.4, 1.6, -1.2], [0, 0, 0, 0], 'ro', markersize=8, label='Предполагаемые корни')
    plt.legend(fontsize=11)
    
    plt.tight_layout()
    plt.show()

def find_intervals():
    """Поиск отрезков с переменой знака"""
    print("\n" + "="*60)
    print("ПОИСК ОТРЕЗКОВ С ПЕРЕМЕНОЙ ЗНАКА")
    print("="*60)
    
    intervals = []
    step = 0.1
    x_start = -2.0
    x_end = 2.0
    
    x = x_start
    while x < x_end:
        try:
            if f(x) * f(x + step) < 0:
                intervals.append((round(x, 2), round(x + step, 2)))
                print(f"Найден отрезок: [{x:.2f}, {x+step:.2f}]")
                print(f"  f({x:.2f}) = {f(x):.6f}")
                print(f"  f({x+step:.2f}) = {f(x+step):.6f}")
        except:
            pass
        x += step
    
    return intervals

def check_newton_conditions(f, df, d2f, a, b, x0):
    """
    Проверка условий сходимости метода Ньютона
    """
    print("ПРОВЕРКА УСЛОВИЙ СХОДИМОСТИ МЕТОДА НЬЮТОНА")
    print("=" * 50)
    
    # 1. Проверка разных знаков на концах отрезка
    fa = f(a)
    fb = f(b)
    print(f"1. f({a}) = {fa:.6f}, f({b}) = {fb:.6f}")
    
    if fa * fb < 0:
        print("✓ f(a)*f(b) < 0 - корень на отрезке есть")
    else:
        print("✗ f(a)*f(b) > 0 - корня на отрезке может не быть")
    
    # 2. Проверка условия |f(x)f''(x)|/(f'(x))² < 1
    print(f"\n2. Проверка условия |f(x)f''(x)|/(f'(x))² < 1:")
    
    # В начальной точке
    for x_check in [x0, a, b]:
        print(f"в точке: {x_check}")
        fx0 = f(x_check)
        dfx0 = df(x_check)
        d2fx0 = d2f(x_check)
        
    
        condition = abs(fx0 * d2fx0) / (dfx0 ** 2)
        print(f"   В {x_check}: {condition:.6f}", end="")
        if condition < 1:
            print(" ✓ < 1")
        else:
            print(" ✗ >= 1")

    print(f"\n3. Знаки производных на отрезке:")
    
    # Первая производная
    df_a = df(a)
    df_b = df(b)
    print(f"   f'({a}) = {df_a:.6f}, f'({b}) = {df_b:.6f}")
    
    if df_a * df_b > 0 and abs(df_a) > eps and abs(df_b) > eps:
        print("   ✓ f'(x) не меняет знак и не равна нулю")
    else:
        print("   ✗ f'(x) меняет знак или равна нулю")
    
    d2f_a = d2f(a)
    d2f_b = d2f(b)
    print(f"   f''({a}) = {d2f_a:.6f}, f''({b}) = {d2f_b:.6f}")
    
    if d2f_a * d2f_b > 0:
        print("   ✓ f''(x) не меняет знак")
    else:
        print("   ✗ f''(x) меняет знак")
    
    print(f"\n4. ИТОГ:")
    
    root_exists = fa * fb < 0
    derivative_nonzero = abs(dfx0) > eps
    derivatives_constant = (df_a * df_b > 0 and d2f_a * d2f_b > 0 and 
                          abs(df_a) > eps and abs(df_b) > eps)
    
    if all([root_exists, derivative_nonzero, derivatives_constant]):
        print("✓ Все условия выполнены - метод Ньютона сойдется")
    else:
        print("✗ Не все условия выполнены - сходимость не гарантирована")

def print_results_table(all_results):
    """Вывод общей таблицы с результатами всех методов для всех отрезков"""
    print("\n\n" + "="*95)
    print("ОБЩАЯ СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ ДЛЯ ВСЕХ ОТРЕЗКОВ")
    print("="*95)
    print(f"{'Отрезок':<15} {'Метод':<20} {'Корень':<15} {'Итерации':<12} {'f(x)':<15}")
    print("-"*95)
    
    for interval, method, root, iters, fx in all_results:
        if iters > 0:
            print(f"{interval:<15} {method:<20} {root:<15.6f} {iters:<12} {fx:<15.10f}")
        else:
            print(f"{interval:<15} {method:<20} {'Не сошёлся':<15} {'-':<12} {'-':<15}")
    print("="*95)

# Проверка условий сходимости
def check_convergence(x0, a, b):
    """Проверка достаточных условий сходимости"""
    print("\n" + "="*60)
    print("ПРОВЕРКА УСЛОВИЙ СХОДИМОСТИ")
    print("="*60)
    
    # Метод Ньютона: f(x0)*f''(x0) > 0 и f'(x) ≠ 0
    print("\n--- Метод Ньютона ---")
    
    check_newton_conditions(f, df, d2f, a, b, x0)
    

# Основная программа
def main():
    eps = 0.0001
    
    print("="*60)
    print("РЕШЕНИЕ НЕЛИНЕЙНОГО УРАВНЕНИЯ (Вариант 41)")
    print("Уравнение: e^(x^2) + 3^(-x) - 5x^2 - 1 = 0")
    print("="*60)
    
    # Поиск отрезков
    intervals = find_intervals()
    
    if not intervals:
        print("Отрезки с переменой знака не найдены!")
        return
    
    i = 0
    all_results = []  # Список для хранения всех результатов
    
    for interval in intervals:
        a, b = interval
        x0 = (a + b) / 2

        print("\n\n" + "="*60)
        print(f"ДЛЯ ОТРЕЗКА: [{a}; {b}]")
        print("\n" + "="*60)
    
        check_convergence(x0, a, b)

        print("\n" + "="*60)
        print("РЕЗУЛЬТАТЫ ВЫЧИСЛЕНИЙ")
        print("="*60)

        results = []
        
        interval_str = f"[{a}; {b}]"
        
        print(f"\nВведите значение lambda для метода простых итераций на отрезке [{a}; {b}]:")
        try:
            lam = float(input("lambda = "))
        except:
            print("Некорректный ввод, используется lambda = 0.1")
            lam = 0.1
        
        # 1. Метод простой итерации
        print("\n1. МЕТОД ПРОСТОЙ ИТЕРАЦИИ")
        try:
            root1, iter1, hist1 = simple_iteration(x0, a, b, lam)
            if len(hist1) > 0:
                print(f"   Корень: x = {root1:.6f}")
                print(f"   Число итераций: {iter1}")
                print(f"   Проверка: f(x) = {f(root1):.10f}")
                all_results.append((interval_str, "Простой итерации", root1, iter1, f(root1)))
        except Exception as e:
            print(f"   Ошибка: {e}")
            all_results.append((interval_str, "Простой итерации", None, None, None))

        # 2. Метод Ньютона
        print("\n2. МЕТОД НЬЮТОНА")
        try:
            root2, iter2, hist2 = newton_method(a, b)
            print(f"   Корень: x = {root2:.6f}")
            print(f"   Число итераций: {iter2}")
            print(f"   Проверка: f(x) = {f(root2):.10f}")
            all_results.append((interval_str, "Ньютона", root2, iter2, f(root2)))
        except Exception as e:
            print(f"   Ошибка: {e}")
            all_results.append((interval_str, "Ньютона", None, None, None))

        # 3. Метод секущих
        print("\n3. МЕТОД СЕКУЩИХ")
        try:
            root3, iter3, hist3 = secant_method(a, b)
            print(f"   Корень: x = {root3:.6f}")
            print(f"   Число итераций: {iter3}")
            print(f"   Проверка: f(x) = {f(root3):.10f}")
            all_results.append((interval_str, "Секущих", root3, iter3, f(root3)))
        except Exception as e:
            print(f"   Ошибка: {e}")
            all_results.append((interval_str, "Секущих", None, None, None))

        # 4. Метод хорд
        print("\n4. МЕТОД ХОРД")
        try:
            root4, iter4, hist4 = chord_method(a, b)
            if root4:
                print(f"   Корень: x = {root4:.6f}")
                print(f"   Число итераций: {iter4}")
                print(f"   Проверка: f(x) = {f(root4):.10f}")
                all_results.append((interval_str, "Хорд", root4, iter4, f(root4)))
            else:
                all_results.append((interval_str, "Хорд", None, None, None))
        except Exception as e:
            print(f"   Ошибка: {e}")
            all_results.append((interval_str, "Хорд", None, None, None))

        # 5. Метод дихотомии
        print("\n5. МЕТОД ДИХОТОМИИ")
        try:
            root5, iter5, hist5 = bisection_method(a, b)
            if root5:
                print(f"   Корень: x = {root5:.6f}")
                print(f"   Число итераций: {iter5}")
                print(f"   Проверка: f(x) = {f(root5):.10f}")
                all_results.append((interval_str, "Дихотомии", root5, iter5, f(root5)))
            else:
                all_results.append((interval_str, "Дихотомии", None, None, None))
        except Exception as e:
            print(f"   Ошибка: {e}")
            all_results.append((interval_str, "Дихотомии", None, None, None))

        i += 1
    
    print_results_table(all_results)

    # Построение графика
    print("\nПостроение графика для определения начального приближения...")
    plot_function()

if __name__ == "__main__":
    main()