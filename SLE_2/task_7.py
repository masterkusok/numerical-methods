import math
import matplotlib.pyplot as plt

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

# 1. Метод простой итерации
def simple_iteration(x0, eps=0.0001, max_iter=1000):
    """Метод простой итерации: x = φ(x)"""
    # Преобразуем: e^(x^2) + 3^(-x) - 5x^2 - 1 = 0
    # => x = ± sqrt((e^(x^2) + 3^(-x) - 1) / 5)
    eps = eps/10
    initial_sign = 1 if x0 >= 0 else -1
    
    def phi(x):
        val = (math.exp(x**2) + 3**(-x) - 1) / 5
        if val < 0:
            return x 
        
        return initial_sign * math.sqrt(val)
    
    # Проверка условия сходимости |φ'(x)| < 1 в окрестности x0
    def phi_derivative_approx(x, h=0.0001):
        return (phi(x + h) - phi(x - h)) / (2 * h)
    
    iterations = []
    x = x0
    iterations.append(x)
    
    # Проверяем условие сходимости
    phi_prime = phi_derivative_approx(x0)
    print(f"   |φ'(x0)| ≈ {abs(phi_prime):.6f}", end="")
    if abs(phi_prime) >= 1:
        print(" - НЕ выполнено условие сходимости!")
    else:
        print(" - условие сходимости выполнено")
        
        for i in range(max_iter):
            x_new = phi(x)
            iterations.append(x_new)
            
            if abs(x_new - x) < eps:
                return x_new, i + 1, iterations
            x = x_new
    
    return x, max_iter, iterations

# 2. Метод Ньютона
def newton_method(x0, eps=0.0001, max_iter=1000):
    """Метод Ньютона (касательных)"""
    iterations = []
    x = x0
    iterations.append(x)
    
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)
        
        if abs(dfx) < 1e-10:
            print("   Производная близка к нулю")
            break
            
        x_new = x - fx / dfx
        iterations.append(x_new)
        
        if abs(x_new - x) < eps:
            return x_new, i + 1, iterations
        x = x_new
    
    return x, max_iter, iterations

# 3. Метод секущих
def secant_method(x0, x1, eps=0.0001, max_iter=1000):
    """Метод секущих"""
    iterations = []
    iterations.append(x0)
    iterations.append(x1)
    
    for i in range(max_iter):
        fx0 = f(x0)
        fx1 = f(x1)
        
        if abs(fx1 - fx0) < 1e-10:
            print("   Деление на ноль в методе секущих")
            break
        
        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        iterations.append(x_new)
        
        if abs(x_new - x1) < eps:
            return x_new, i + 1, iterations
        
        x0, x1 = x1, x_new
    
    return x1, max_iter, iterations
# 4. Метод хорд
def chord_method(a, b, eps=0.0001, max_iter=1000):
    """Метод хорд (метод ложного положения)"""
    iterations = []
    
    fa = f(a)
    fb = f(b)
    
    print(f"   f({a}) = {fa:.6f}, f({b}) = {fb:.6f}")
    
    if fa * fb > 0:
        print("   ОШИБКА: f(a)*f(b) > 0 - корня на отрезке нет!")
        return None, 0, []
    
    # Определяем неподвижный конец (где знак f(x) совпадает со знаком f''(x))
    d2fa = d2f(a)
    d2fb = d2f(b)
    
    print(f"   f''({a}) = {d2fa:.6f}, f''({b}) = {d2fb:.6f}")
    
    if fa * d2fa > 0:
        print(f"   Неподвижный конец: a = {a}")
        x = b
        iterations.append(x)
        for i in range(max_iter):
            fx = f(x)
            if abs(fx - fa) < 1e-10:
                break
            x_new = x - fx * (x - a) / (fx - fa)
            iterations.append(x_new)
            
            if abs(x_new - x) < eps:
                return x_new, i + 1, iterations
            x = x_new
    else:
        print(f"   Неподвижный конец: b = {b}")
        x = a
        iterations.append(x)
        for i in range(max_iter):
            fx = f(x)
            if abs(fx - fb) < 1e-10:
                break
            x_new = x - fx * (x - b) / (fx - fb)
            iterations.append(x_new)
            
            if abs(x_new - x) < eps:
                return x_new, i + 1, iterations
            x = x_new
    
    return x, max_iter, iterations

# 5. Метод дихотомии (деления пополам)
def bisection_method(a, b, eps=0.0001, max_iter=1000):
    """Метод дихотомии"""
    iterations = []
    
    fa = f(a)
    fb = f(b)
    
    print(f"   f({a}) = {fa:.6f}, f({b}) = {fb:.6f}")
    
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

# Графическое определение начального приближения
def plot_function():
    """Построение графика функции"""
    x_vals = [i * 0.01 for i in range(-150, 151)]
    y_vals = []
    
    for x in x_vals:
        try:
            y = f(x)
            if abs(y) < 50:  # Ограничиваем для читаемости графика
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
    
    # Отмечаем предполагаемые корни
    plt.plot([-0.85, 0.4], [0, 0], 'ro', markersize=8, label='Предполагаемые корни')
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
    x_start = -1.0
    x_end = 1.0
    
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
    fx0 = f(x0)
    dfx0 = df(x0)
    d2fx0 = d2f(x0)
    
    if abs(dfx0) > 1e-10:
        condition1 = abs(fx0 * d2fx0) / (dfx0 ** 2)
        print(f"   В x0 = {x0}: {condition1:.6f}", end="")
        if condition1 < 1:
            print(" ✓ < 1")
        else:
            print(" ✗ >= 1")
    else:
        print("   ✗ f'(x0) слишком близко к нулю")
    
    # На всем отрезке
    max_condition = 0
    for i in range(10):  # Проверяем в 10 точках
        x = a + i * (b - a) / 9
        fx = f(x)
        dfx = df(x)
        d2fx = d2f(x)
        
        if abs(dfx) > 1e-10:
            cond = abs(fx * d2fx) / (dfx ** 2)
            if cond > max_condition:
                max_condition = cond
    
    print(f"   Максимум на отрезке: {max_condition:.6f}", end="")
    if max_condition < 1:
        print(" ✓ < 1")
    else:
        print(" ✗ >= 1")
    
    # 3. Проверка постоянства знаков производных
    print(f"\n3. Знаки производных на отрезке:")
    
    # Первая производная
    df_a = df(a)
    df_b = df(b)
    print(f"   f'({a}) = {df_a:.6f}, f'({b}) = {df_b:.6f}")
    
    if df_a * df_b > 0 and abs(df_a) > 1e-10 and abs(df_b) > 1e-10:
        print("   ✓ f'(x) не меняет знак и не равна нулю")
    else:
        print("   ✗ f'(x) меняет знак или равна нулю")
    
    # Вторая производная
    d2f_a = d2f(a)
    d2f_b = d2f(b)
    print(f"   f''({a}) = {d2f_a:.6f}, f''({b}) = {d2f_b:.6f}")
    
    if d2f_a * d2f_b > 0:
        print("   ✓ f''(x) не меняет знак")
    else:
        print("   ✗ f''(x) меняет знак")
    
    # 4. Итоговый вывод
    print(f"\n4. ИТОГ:")
    
    # Проверяем все условия
    root_exists = fa * fb < 0
    derivative_nonzero = abs(dfx0) > 1e-10
    condition_satisfied = max_condition < 1
    derivatives_constant = (df_a * df_b > 0 and d2f_a * d2f_b > 0 and 
                          abs(df_a) > 1e-10 and abs(df_b) > 1e-10)
    
    if all([root_exists, derivative_nonzero, condition_satisfied, derivatives_constant]):
        print("✓ Все условия выполнены - метод Ньютона сойдется")
        print("✓ Так как сходится метод Ньютона - то сходятся и методы секущих и хорд.")
    else:
        print("✗ Не все условия выполнены - сходимость не гарантирована")

# Проверка условий сходимости
def check_convergence(x0, a, b):
    """Проверка достаточных условий сходимости"""
    print("\n" + "="*60)
    print("ПРОВЕРКА УСЛОВИЙ СХОДИМОСТИ")
    print("="*60)
    
    try:
        fx0 = f(x0)
        dfx0 = df(x0)
        d2fx0 = d2f(x0)
    except:
        print("Ошибка вычисления в точке x0")
        return
    
    print(f"\nНачальное приближение: x0 = {x0}")
    print(f"Отрезок локализации: [{a}, {b}]")
    print(f"f(x0) = {fx0:.6f}")
    print(f"f'(x0) = {dfx0:.6f}")
    print(f"f''(x0) = {d2fx0:.6f}")
    
    # Метод Ньютона: f(x0)*f''(x0) > 0 и f'(x) ≠ 0
    print("\n--- Метод Ньютона ---")
    
    check_newton_conditions(f, df, d2f, a, b, x0)
    

# Основная программа
def main():
    eps = 0.0001
    
    print("="*60)
    print("РЕШЕНИЕ НЕЛИНЕЙНОГО УРАВНЕНИЯ (Вариант 41)")
    print("Уравнение: e^(x^2) + 3^(-x) - 5x^2 - 1 = 0")
    print(f"Точность: ε = {eps}")
    print("="*60)
    
    # Построение графика
    print("\nПостроение графика для определения начального приближения...")
    plot_function()
    
    # Поиск отрезков
    intervals = find_intervals()
    
    if not intervals:
        print("Отрезки с переменой знака не найдены!")
        return
    
    # Выбираем положительный корень
    print("\n" + "="*60)
    print("РАБОТАЕМ С ПОЛОЖИТЕЛЬНЫМ КОРНЕМ")
    print("="*60)
    
    
    
    for interval in intervals:
        a, b = interval
        x0 = (a + b) / 2
        print("\n\n\n\n\n\n" + str(x0) + "\n\n\n\n\n\n\n\n")

        print("\n\n" + "="*60)
        print(f"ДЛЯ ОТРЕЗКА: [{a}; {b}]")
        print("\n" + "="*60)
    
        check_convergence(x0, a, b)

        print("\n" + "="*60)
        print("РЕЗУЛЬТАТЫ ВЫЧИСЛЕНИЙ")
        print("="*60)

        results = []

        # 1. Метод простой итерации
        print("\n1. МЕТОД ПРОСТОЙ ИТЕРАЦИИ")
        try:
            root1, iter1, hist1 = simple_iteration(x0, eps)
            print(f"   Корень: x = {root1:.6f}")
            print(f"   Число итераций: {iter1}")
            print(f"   Проверка: f(x) = {f(root1):.10f}")
            results.append(("Простой итерации", root1, iter1, f(root1), hist1))
        except Exception as e:
            print(f"   Ошибка: {e}")
            results.append(("Простой итерации", None, None, None, []))

        # 2. Метод Ньютона
        print("\n2. МЕТОД НЬЮТОНА")
        try:
            root2, iter2, hist2 = newton_method(x0, eps)
            print(f"   Корень: x = {root2:.6f}")
            print(f"   Число итераций: {iter2}")
            print(f"   Проверка: f(x) = {f(root2):.10f}")
            results.append(("Ньютона", root2, iter2, f(root2), hist2))
        except Exception as e:
            print(f"   Ошибка: {e}")
            results.append(("Ньютона", None, None, None, []))

        # 3. Метод секущих
        print("\n3. МЕТОД СЕКУЩИХ")
        try:
            root3, iter3, hist3 = secant_method(a, b, eps)
            print(f"   Корень: x = {root3:.6f}")
            print(f"   Число итераций: {iter3}")
            print(f"   Проверка: f(x) = {f(root3):.10f}")
            results.append(("Секущих", root3, iter3, f(root3), hist3))
        except Exception as e:
            print(f"   Ошибка: {e}")
            results.append(("Секущих", None, None, None, []))


        # 4. Метод хорд
        print("\n4. МЕТОД ХОРД")
        try:
            root4, iter4, hist4 = chord_method(a, b, eps)
            if root4:
                print(f"   Корень: x = {root4:.6f}")
                print(f"   Число итераций: {iter4}")
                print(f"   Проверка: f(x) = {f(root4):.10f}")
                results.append(("Хорд", root4, iter4, f(root4), hist4))
            else:
                results.append(("Хорд", None, None, None, []))
        except Exception as e:
            print(f"   Ошибка: {e}")
            results.append(("Хорд", None, None, None, []))

        # 5. Метод дихотомии
        print("\n5. МЕТОД ДИХОТОМИИ")
        try:
            root5, iter5, hist5 = bisection_method(a, b, eps)
            if root5:
                print(f"   Корень: x = {root5:.6f}")
                print(f"   Число итераций: {iter5}")
                print(f"   Проверка: f(x) = {f(root5):.10f}")
                results.append(("Дихотомии", root5, iter5, f(root5), hist5))
            else:
                results.append(("Дихотомии", None, None, None, []))
        except Exception as e:
            print(f"   Ошибка: {e}")
            results.append(("Дихотомии", None, None, None, []))

if __name__ == "__main__":
    main()